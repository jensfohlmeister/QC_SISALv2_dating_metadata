
using Pkg
#Pkg.add("Conda")
#Pkg.build("PyCall")
#Pkg.add("PyCall")
Pkg.add("Plots")
Pkg.add("Statistics")
using Plots
using Statistics


# List of input variables:
# UU = measured 234U/238U activity ratio
# dUU = measured 2sigma error of UU
# Th0U = measured 230Th/238U activity ratio
# dTh0U = measured 2 sigma error of Th0U
# Th2U = measured 232Th/238U activity ratio
# dTh2U = measured 2 sigma error of Th2U
# Th0Th2ini = assumed initial 230Th/232Th activity ratio
# dTh0Th2ini = assumed 2 sigma error of Th0Th2ini

# List of resulting variables
# t_uncorr = uncorrected age [a before measurement year]
# dt_uncorr = error of uncorrected age [a]
# t_corr = detrital Th corrected age [a before measurement year]
# dt_corr = error of detrital Th corrected age [a]

#Corrected 230Th ages assume the initial 230Th/232Th atomic ratio of
#4.4  ±2.2 x10-6. This equals 230Th/232Th activity ratio of 0.818+/-0.409.
#Those are the values for a material at secular
#equilibrium, with the bulk earth 232Th/238U value of 3.8.
#The errors are arbitrarily assumed to be 50%.


###############################################################################
function age_determination(UU1, ThU1, Decay0, Decay4)

    t1 = 0
    t2 = 1000000

    x1 = (1 - exp(-t1 * Decay0)) + (UU1 - 1) * Decay0 / (Decay0 - Decay4) *
            (1-exp(-t1 * (Decay0 - Decay4))) - ThU1
    x2 = (1 - exp(-t2 * Decay0)) + (UU1 - 1) * Decay0 / (Decay0 - Decay4) *
            (1-exp(-t2 * (Decay0 - Decay4))) -ThU1

    if x1*x2 >= 0
        t = "out of range"

    else
        while abs(t2-t1) > 0.01
            t = (t1 + t2) / 2
            res = (1 - exp(-t * Decay0)) + (UU1 - 1) * Decay0 / (Decay0 - Decay4) *
                (1-exp(-t * (Decay0 - Decay4))) - ThU1
            if res > 0
                t2 = t
            else
                t1 = t
            end
        end
    end

    return t
end

###############################################################################
function age_err_determination(UU1,dUU1,ThU1,dThU1, Decay0, Decay4)

    t = age_determination(UU1, ThU1, Decay0, Decay4)

    nrMC = 5000
    t1 = zeros(nrMC)
    UUr = randn(nrMC)*dUU1/2 .+ UU1
    ThUr = randn(nrMC)*dThU1/2 .+ ThU1

    for i = 1:nrMC
        tdummy = age_determination(UUr[i],ThUr[i], Decay0, Decay4)

        if tdummy == "out of range"
            i=i-1
        else
            t1[i] = tdummy
        end
    end

    return t, 2*std(t1)
end

###############################################################################
function detritus_corr(UU, Th0U, Th2U, Th0Th2ini, Decay0, Decay4)

    t1 = 0
    t2 = 1000000

    x1 = (1 - exp(-t1 * Decay0)) + Th2U * Th0Th2ini * exp(-t1 * Decay0) +
            (UU - 1) * Decay0 / (Decay0 - Decay4) *
            (1-exp(-t1 * (Decay0 - Decay4))) - Th0U
    x2 = (1 - exp(-t2 * Decay0)) + Th2U * Th0Th2ini * exp(-t2 * Decay0) +
            (UU - 1) * Decay0 / (Decay0 - Decay4) *
            (1-exp(-t2 * (Decay0 - Decay4))) - Th0U

    if x1*x2 >= 0
        t = "out of range"

    else
        while abs(t2-t1) > 0.01
            t = (t1 + t2) / 2
            res = (1 - exp(-t * Decay0)) + Th2U * Th0Th2ini * exp(-t * Decay0) +
                (UU - 1) * Decay0 / (Decay0 - Decay4) *
                (1-exp(-t * (Decay0 - Decay4))) - Th0U
            if res > 0
                t2 = t
            else
                t1 = t
            end
        end
    end

    return t
end

###############################################################################
function detritus_corr_err(UU,dUU,Th0U,dTh0U,Th2U,dTh2U,Th0Th2ini,dTh0Th2ini, Decay0, Decay4)

    t = detritus_corr(UU, Th0U, Th2U, Th0Th2ini, Decay0, Decay4)

    nrMC = 5000
    t1 = zeros(nrMC)
    UUr = randn(nrMC)*dUU/2 .+ UU
    Th0Ur = randn(nrMC)*dTh0U/2 .+ Th0U
    Th2Ur = randn(nrMC)*dTh2U/2 .+ Th2U
    Th0Th2inir = randn(nrMC)*dTh0Th2ini/2 .+ Th0Th2ini

    for i = 1:nrMC

        tdummy = detritus_corr(UUr[i], Th0Ur[i], Th2Ur[i], Th0Th2inir[i], Decay0, Decay4)

        if tdummy == "out of range"
            i=i-1
        else
            t1[i] = tdummy
        end
    end

    return t, 2*std(t1)
end

###############################################################################


################################################################################
# preparation of variables and errors for dating procedure
################################################################################
function prepare1(var,var_err,rel_err)
    # var: variable; var_err: variables error; rel_err: alternatively use relative error
    if ismissing(var_err)
        var_err1 = var*rel_err  #assume 1% error
    else
        var_err1 = var_err  #use given error
    end
    return (var, var_err1)
end
function prepare2(var1,var1_err,rel_err1,var2,var2_err,rel_err2)
    # var: variable; var_err: variables error; rel_err: alternatively use relative error
    var = var1/var2
    if ismissing(var1_err) & ismissing(var2_err)
        var_err = sqrt((var1*rel_err1)^2 + (var2*rel_err2)^2)
    elseif ismissing(var1_err)
        var_err = sqrt((var1*rel_err1)^2 + (var1*var2_err/var2^2)^2)
    elseif ismissing(var2_err)
        var_err = sqrt((var1_err/var2)^2 + (var2*rel_err2)^2)
    else
        var_err = sqrt((var1_err/var2)^2 + (var1*var2_err/var2^2)^2)  #use given error
    end
    return (var, var_err)
end
################################################################################
