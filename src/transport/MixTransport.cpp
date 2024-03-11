/**
 *  @file MixTransport.cpp
 *  Mixture-averaged transport properties for ideal gas mixtures.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/transport/MixTransport.h"
#include "cantera/base/stringUtils.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/base/utilities.h"

namespace Cantera
{

void MixTransport::init(ThermoPhase* thermo, int mode, int log_level)
{
    GasTransport::init(thermo, mode, log_level);
    m_cond.resize(m_nsp);
    // precompute and store log(epsilon_ij/k_B)
    m_log_eps_k.resize(m_nsp, m_nsp);
    for (size_t i = 0; i < m_nsp; i++) {
        for (size_t j = i; j < m_nsp; j++) {
            m_log_eps_k(i,j) = log(m_epsilon(i,j)/Boltzmann);
            m_log_eps_k(j,i) = m_log_eps_k(i,j);
        }
    }
    m_astar.resize(m_nsp, m_nsp);
    m_bstar.resize(m_nsp, m_nsp);
    m_cstar.resize(m_nsp, m_nsp);

}

void MixTransport::getMobilities(double* const mobil)
{
    getMixDiffCoeffs(m_spwork.data());
    double c1 = ElectronCharge / (Boltzmann * m_temp);
    for (size_t k = 0; k < m_nsp; k++) {
        mobil[k] = c1 * m_spwork[k];
    }
}

double MixTransport::thermalConductivity()
{
    update_T();
    update_C();
    if (!m_spcond_ok) {
        updateCond_T();
    }
    if (!m_condmix_ok) {
        double sum1 = 0.0, sum2 = 0.0;
        for (size_t k = 0; k < m_nsp; k++) {
            sum1 += m_molefracs[k] * m_cond[k];
            sum2 += m_molefracs[k] / m_cond[k];
        }
        m_lambda = 0.5*(sum1 + 1.0/sum2);
        m_condmix_ok = true;
    }
    return m_lambda;
}

void MixTransport::getThermalDiffCoeffs(double* const dt)
{

    //
    // Update the temperature-dependent quantities
    //
    double print_counter = 0.0;
    update_T(); // Update the temperature-dependent quantities
    update_C(); // Update the mass fractions
    updateSpeciesViscosities(); // Update the species viscosities
    updateDiff_T();  // Update the binary diffusive coefficients
    getMixDiffCoeffs(m_spwork.data());  // Update the mixture-averaged diffusion coefficients
    
    // evaluate polynomial fits for A*, B*, C*
    for (size_t i = 0; i < m_nsp; i++) {
        for (size_t j = i; j < m_nsp; j++) {
            double z = m_star_poly_uses_actualT[i][j] == 1 ? m_logt : m_logt - m_log_eps_k(i,j);
            int ipoly = m_poly[i][j];
            if (m_mode == CK_Mode) {
                m_astar(i,j) = poly6(z, m_astar_poly[ipoly].data());
                m_bstar(i,j) = poly6(z, m_bstar_poly[ipoly].data());
                m_cstar(i,j) = poly6(z, m_cstar_poly[ipoly].data());
            } else {
                m_astar(i,j) = poly8(z, m_astar_poly[ipoly].data());
                m_bstar(i,j) = poly8(z, m_bstar_poly[ipoly].data());
                m_cstar(i,j) = poly8(z, m_cstar_poly[ipoly].data());
            }
            m_astar(j,i) = m_astar(i,j);
            m_bstar(j,i) = m_bstar(i,j);
            m_cstar(j,i) = m_cstar(i,j);
        }
    }

    const double* Y_i = m_thermo->massFractions();
    double p = m_thermo->pressure();
    double mean_w = m_thermo->meanMolecularWeight();
    double X_i[m_nsp]; // Not sure if m_molefracs is already defined
    double rho = m_thermo->density();
    double S;  // Mass conservation term
    S = 0.0;       
    double mu;
    mu = 0.0;
    for (size_t k = 0; k < m_nsp; k++) {
        X_i[k] = Y_i[k] * mean_w / m_mw[k];
        mu += X_i[k] * m_visc[k];
    }
       
    double kT_i[m_nsp];  // k^T_i
    double a_i[m_nsp];  // a_k
    double lambda_mon_i[m_nsp];  // lambda^mon_i
    double Delta_i[m_nsp];  // Delta_i
    double Phi_ij[m_nsp][m_nsp];  // Phi_ij

    for (size_t k = 0; k < m_nsp; k++) {
        for (size_t j = 0; j < m_nsp; j++) {
            Phi_ij[k][j] = 0.0;
            Phi_ij[k][j] = (1.065/(2*sqrt(2))) * pow(1 + pow(m_visc[k]/m_visc[j], 0.5) * pow(m_mw[j]/m_mw[k], 0.25), 2);
            Phi_ij[k][j] = Phi_ij[k][j] / pow(1 + m_mw[k]/m_mw[j], 0.5);
        }
    }

    for (size_t k = 0; k < m_nsp; k++) {
        Delta_i[k] = X_i[k];
        for (size_t j = 0; j < m_nsp; j++) {
            Delta_i[k] += X_i[j] * Phi_ij[k][j];
        }
        Delta_i[k] -= X_i[k] * Phi_ij[k][k];  // Remove self-contribution
    }

    for (size_t k = 0; k < m_nsp; k++) {
        lambda_mon_i[k] = (15.0/4.0) * (GasConstant * m_visc[k] / m_mw[k]);
        a_i[k] = lambda_mon_i[k] * X_i[k] / Delta_i[k];
    }

    // Compute thermal diffusion coefficients D^T_k
    for (size_t k = 0; k < m_nsp; k++) {
        kT_i[k] = 0.0;
        for (size_t j = 0; j < m_nsp; j++) {
            kT_i[k] += (1.2 * m_cstar(k,j) - 1) * (Y_i[k]*a_i[j] - Y_i[j]*a_i[k]) / ((m_bdiff(k,j)/p)*(m_mw[k] + m_mw[j]));
        }
        kT_i[k] = kT_i[k] * mean_w * mean_w / GasConstant / rho;
    }

    // Compute thermal diffusion coefficients D^T_k
    for (size_t k = 0; k < m_nsp; k++) {
        dt[k] = rho * m_mw[k] * kT_i[k] * m_spwork[k] / mean_w;
        S += dt[k];
    }

}

void MixTransport::getSpeciesFluxes(size_t ndim, const double* const grad_T,
                                    size_t ldx, const double* const grad_X,
                                    size_t ldf, double* const fluxes)
{
    update_T();
    update_C();
    getMixDiffCoeffs(m_spwork.data());
    const vector<double>& mw = m_thermo->molecularWeights();
    const double* y = m_thermo->massFractions();
    double rhon = m_thermo->molarDensity();
    vector<double> sum(ndim,0.0);
    for (size_t n = 0; n < ndim; n++) {
        for (size_t k = 0; k < m_nsp; k++) {
            fluxes[n*ldf + k] = -rhon * mw[k] * m_spwork[k] * grad_X[n*ldx + k];
            sum[n] += fluxes[n*ldf + k];
        }
    }
    // add correction flux to enforce sum to zero
    for (size_t n = 0; n < ndim; n++) {
        for (size_t k = 0; k < m_nsp; k++) {
            fluxes[n*ldf + k] -= y[k]*sum[n];
        }
    }
}

void MixTransport::update_T()
{
    double t = m_thermo->temperature();
    if (t == m_temp && m_nsp == m_thermo->nSpecies()) {
        return;
    }
    if (t < 0.0) {
        throw CanteraError("MixTransport::update_T",
                           "negative temperature {}", t);
    }
    GasTransport::update_T();
    // temperature has changed, so polynomial fits will need to be redone.
    m_spcond_ok = false;
    m_bindiff_ok = false;
    m_condmix_ok = false;
}

void MixTransport::update_C()
{
    // signal that concentration-dependent quantities will need to be recomputed
    // before use, and update the local mole fractions.
    m_visc_ok = false;
    m_condmix_ok = false;
    m_thermo->getMoleFractions(m_molefracs.data());

    // add an offset to avoid a pure species condition
    for (size_t k = 0; k < m_nsp; k++) {
        m_molefracs[k] = std::max(Tiny, m_molefracs[k]);
    }
}

void MixTransport::updateCond_T()
{
    if (m_mode == CK_Mode) {
        for (size_t k = 0; k < m_nsp; k++) {
            m_cond[k] = exp(dot4(m_polytempvec, m_condcoeffs[k]));
        }
    } else {
        for (size_t k = 0; k < m_nsp; k++) {
            m_cond[k] = m_sqrt_t * dot5(m_polytempvec, m_condcoeffs[k]);
        }
    }
    m_spcond_ok = true;
    m_condmix_ok = false;
}

}
