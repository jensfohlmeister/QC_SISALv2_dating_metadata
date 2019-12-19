
# Comment blocks: 'ctrl' + '/'

using Pkg

Pkg.add("CSV")
Pkg.add("DataFrames")
Pkg.add("DataFramesMeta")


using CSV
using DataFrames, DataFramesMeta, Logging

include("U-ThAgeCalc_function.jl")


################################################################################
# define global variables
################################################################################
global l234new, l234old, l230old, l230new, l230, l234

U234_Thalf_new = 245620 # Cheng et al 2013, EPSL
Th230_Thalf_new = 75584 # Cheng et al 2013, EPSL
U234_Thalf_old = 245250 # Cheng et al 2000, Chem Geol
Th230_Thalf_old = 75690 # Cheng et al 2000, Chem Geol

l234new = log(2)/U234_Thalf_new
l234old = log(2)/U234_Thalf_old
l234veryold = 2.835*10^(-6)     # Edwards et al. 1987, EPSL
l230veryold = 9.195*10^(-6)     # Edwards et al. 1987, EPSL
l230old = log(2)/Th230_Thalf_old
l230new = log(2)/Th230_Thalf_new

l232 = 0.000000000049475 #232Th decay constant; calculated with l values from CRC Handbook
l238 = 0.000000000155125 #238U decay constant; Jaffey et al., 1971
################################################################################

################################################################################
# reading data files
################################################################################

cd(@__DIR__) # to change to path of this *.jl file
cd("SISALv2_csv") # to change the path relative to where it was before
entity = CSV.read("entity.csv",missingstring = "NULL")
gap = CSV.read("gap.csv",types=(Int64,String)) # reads CSV file relative to where the julia is running
            # Surprisingly this is not the same, where your *.jl file is located
composite_link_entity = CSV.read("composite_link_entity.csv")
d13C = CSV.read("d13C.csv",missingstring = "NULL")
d18O = CSV.read("d18O.csv",missingstring = "NULL")
dating_lamina = CSV.read("dating_lamina.csv")
dating = CSV.read("dating.csv",missingstring = "NULL")#treating all NULL as missing values
names!(dating,Symbol.([:dating_id,
    :entity_id,:date_type,:depth_dating,:dating_thickness,:lab_num,:material_dated,
    :min_weight,:max_weight,:uncorr_age,:uncorr_age_uncert_pos,:uncorr_age_uncert_neg,
    :C14_correction,:calib_used,:date_used,:c238U_content,:c238U_uncertainty,
    :c232Th_content,:c232Th_uncertainty,:c230Th_content,:c230Th_uncertainty,
    :a230Th_232Th_ratio,:a230Th_232Th_ratio_uncertainty,:a230Th_238U_activity,
    :a230Th_238U_activity_uncertainty,:a234U_238U_activity,:a234U_238U_activity_uncertainty,
    :ini_230Th_232Th_ratio,:ini_230Th_232Th_ratio_uncertainty,:decay_constant,
    :corr_age,:corr_age_uncert_pos,:corr_age_uncert_neg,:date_used_COPRA,:date_used_linear,
    :date_used_linear_regress,:date_used_Bchron,:date_used_Bacon]))
    # it is necessary to rename those lines with a number on first position
entity_link_reference = CSV.read("entity_link_reference.csv")
hiatus = CSV.read("hiatus.csv")
notes = CSV.read("notes.csv")
original_chronology = CSV.read("original_chronology.csv")
reference = CSV.read("reference.csv")
sample = CSV.read("sample.csv")
sisal_chronology = CSV.read("sisal_chronology.csv")
site = CSV.read("site.csv")
cd(@__DIR__)
################################################################################


################################################################################
#### find all U/Th ages (neglect 14C and events, neglect dirty subsamples)
####                   ('missing' in decay-constant, 'missing' in corr_age)
#### find all U/Th ages in [10000,120000] --> datingDO
################################################################################
datingDO = dropmissing(dating, :decay_constant)#,:corr_age)
#datingDO = datingDO[10000 .<= datingDO.corr_age .<= 120000,:]

### all stalagmites in dating.csv
stalDO = sort(unique(datingDO,:entity_id),:entity_id)

### all stalagmites with DOs in entity.csv
#entity1=stalDO[!,:entity_id]
entityDO = entity[findall(in(stalDO[!,:entity_id]),entity.entity_id),[:site_id,:entity_id,:entity_name,:contact,:data_DOI_URL]]
println(size(entityDO,1)," stalagmites avaialable in SISAL covering the requested period.")

### all cave sites with stalagmites with DOs
siteDO = site[findall(in(entityDO[!,:site_id]),site.site_id),:]
println("Those stalagmites are from $(size(siteDO,1)) caves around the world.")

# output for table in paper or for a world map within the paper about
# stalagmite location/spatial coverage realised in another software
siteDO1 = site[findall(in(entityDO[!,:site_id]),site.site_id),[:site_id, :site_name, :latitude, :longitude, :elevation]]
#CSV.write("siteDO1.csv", siteDO1)
################################################################################

referenceDO = entity_link_reference[findall(in(entityDO[!,:entity_id]),entity_link_reference.entity_id),:]
referenceDO1 = reference[findall(in(referenceDO[!,:ref_id]),reference.ref_id),:]
#datingDO1 = dropmissing(datingDO, :a234U_238U_activity)

df = DataFrame() # output 'df' required by Laia
dummy = DataFrame(site_name = String[], latitude = Float64[], longitude = Float64[] )
df.site_id = entityDO.site_id
for i = 1:size(df.site_id,1)
    for j = 1:size(siteDO.site_id,1)
        if df.site_id[i]==siteDO.site_id[j]
            push!(dummy,(siteDO.site_name[j],siteDO.latitude[j],siteDO.longitude[j]))
        end
    end
end

df.site_name = dummy.site_name
df.entity_id = entityDO.entity_id
df.entity_name = entityDO.entity_name
df.longitude = dummy.longitude
df.latitude = dummy.latitude

dummy1 = DataFrame()
dummy1.ref_id = 1:size(entityDO.entity_id,1)
for i = 1:size(df.site_id,1)
    for j = 1:size(referenceDO.entity_id,1)
        if df.entity_id[i]==referenceDO.entity_id[j]
            dummy1.ref_id[i]=referenceDO.ref_id[j] #another way of filling the
                    #column for a DataFrame
        end
    end
end

dummy2 = DataFrame(citation = String[], publication_DOI = String[])
for i = 1:size(dummy1.ref_id,1)
    for j = 1:size(referenceDO1.ref_id,1)
        if dummy1.ref_id[i]==referenceDO1.ref_id[j]
            push!(dummy2,(referenceDO1.citation[j],referenceDO1.publication_DOI[j]))
        end
    end
end

df.publication_DOI = dummy2.publication_DOI
df.citation = dummy2.citation
df.data_DOI_URL = entityDO.data_DOI_URL
df.contact = entityDO.contact

################################################################################
### starting rough quality control of U-Th dating sheet
################################################################################


c238U = datingDO.c238U_content
c232Th = datingDO.c232Th_content
c230Th = datingDO.c230Th_content
a230Th232Th = datingDO.a230Th_232Th_ratio
a230Th238U = datingDO.a230Th_238U_activity
a234U238U = datingDO.a234U_238U_activity
a230Th232Thini = datingDO.ini_230Th_232Th_ratio
dating_idx = datingDO.dating_id
entity_idx = datingDO.entity_id

# ################################################################################
# ### following the code on "https://docs.julialang.org/en/v1/stdlib/Logging/#Writing-log-events-to-a-file-1"
# ### on "Logging" --> "Writing log events to a file"
# ### all necessary lines are marked by "#log" - NOT USED HERE
# ################################################################################
#
# io = open("log_234U238U.txt", "w+") #log
# logger = SimpleLogger(io) #log
#
# ### check for column a234U238U
#
# entity_idx_no234U238U = []#zeros(size(c238U,1))
# entity_idx_d234U = []
# entity_idx_possible_d234U = []
#
# for i = 1:size(c238U,1)
#     if ismissing(a234U238U[i])==true
#         push!(entity_idx_no234U238U, entity_idx[i])
#         with_logger(logger) do #log
#             @info("no 234U/238U activity is given for dating_id $(dating_idx[i])
#             Please check for all 234U/238 data of this entity ($(entity_idx[i]))")
#         end
#     elseif a234U238U[i]<=0
#         push!(entity_idx_d234U, entity_idx[i])
#         with_logger(logger) do #log
#             @info("d234U is inserted instead of activity ratio for dating_id $(dating_idx[i])
#             Please check for all 234U/238 data of this entity ($(entity_idx[i]))")
#         end
#     elseif a234U238U[i]>10
#         push!(entity_idx_possible_d234U, entity_idx[i])
#         with_logger(logger) do #log
#             @info("most likely d234U is inserted instead of activity ratio for dating_id $(dating_idx[i])
#             Please check for all 234U/238 data of this entity ($(entity_idx[i]))")
#         end
#     end
# end
#
# unique!(entity_idx_no234U238U)
# unique!(entity_idx_d234U)
# unique!(entity_idx_possible_d234U)
#
# flush(io) #log
# close(io) #log


################################################################################
### check wrong values
################################################################################
function checkUTh(checkvar, entity_idx; limit1=missing, limit2=missing)
    ### checks if information is missing and stores entity_id
    checkvar_out = []
    checkvar_out1 = []
    checkvar_out2 = []

    for i = 1:size(checkvar,1)
        if ismissing(checkvar[i])==true
            push!(checkvar_out, entity_idx[i])
        elseif ismissing(limit1)==false && checkvar[i]<=limit1
            push!(checkvar_out1, entity_idx[i])
        elseif ismissing(limit2)==false && checkvar[i]>=limit2
            push!(checkvar_out2, entity_idx[i])
        end
    end

    unique!(checkvar_out)
    unique!(checkvar_out1)
    unique!(checkvar_out2)

    if ismissing(limit1)==false && ismissing(limit2)==false
        return checkvar_out, checkvar_out1, checkvar_out2
    elseif ismissing(limit1)==false
        return checkvar_out,checkvar_out1
    elseif ismissing(limit2)==false
        return checkvar_out,checkvar_out2
    else
        return checkvar_out
    end
end

################################################################################

entity_idx_no238U = checkUTh(c238U,entity_idx)
entity_idx_no232Th = checkUTh(c232Th,entity_idx)
entity_idx_no230Th = checkUTh(c230Th,entity_idx)
entity_idx_no234U238U = checkUTh(a234U238U,entity_idx,limit1 = 0, limit2 = 10)
entity_idx_no230Th238U = checkUTh(a230Th238U,entity_idx)
entity_idx_no230Th232Th = checkUTh(a230Th232Th,entity_idx, limit1 = 0.1)
entity_idx_no230Th232Thini = checkUTh(a230Th232Thini,entity_idx, limit1 = 0.1)



################################################################################
### add notes to each stalagmite if data are missing
################################################################################
function check_entity(df,entity_idx,label)
    dummy3 = DataFrame(missing_x = String[])
    for i = 1:size(df.entity_id,1)
        ii = findall(x -> x==df.entity_id[i],entity_idx)
        if isempty(ii)
            push!(dummy3.missing_x, "$(label) available")
        else
            push!(dummy3.missing_x, "no $(label)")
        end
    end
    return dummy3
end

dummy3 = check_entity(df, entity_idx_no238U, "238U")
df.c238U = dummy3.missing_x
dummy3 = check_entity(df, entity_idx_no232Th, "232Th")
df.c232Th = dummy3.missing_x
dummy3 = check_entity(df, entity_idx_no230Th, "230Th")
df.c230Th = dummy3.missing_x
dummy3 = check_entity(df, entity_idx_no230Th232Th[1], "a230Th/232Th")
dummy4 = check_entity(df, entity_idx_no230Th232Th[2], "no a230Th/232Th but atomic ratio(?)")
for i = 1:size(entity_idx_no230Th232Th[2],1)
    dummy3.missing_x[findall(x->x==entity_idx_no230Th232Th[2][i],df.entity_id)] =
        dummy4.missing_x[findall(x->x==entity_idx_no230Th232Th[2][i],df.entity_id)]
end
df.a230Th232Th = dummy3.missing_x
dummy3 = check_entity(df, entity_idx_no230Th238U, "230Th/238U")
df.a230Th238U = dummy3.missing_x
dummy3 = check_entity(df, entity_idx_no234U238U[1], "234U/238U")
dummy4 = check_entity(df, entity_idx_no234U238U[2], "234U/238U, but d234U")
for i = 1:size(entity_idx_no234U238U[2],1)
    dummy3.missing_x[findall(x->x==entity_idx_no234U238U[2][i],df.entity_id)] =
        dummy4.missing_x[findall(x->x==entity_idx_no234U238U[2][i],df.entity_id)]
end
dummy4 = check_entity(df, entity_idx_no234U238U[3], "a234U/238U, but d234U(?)")
for i = 1:size(entity_idx_no234U238U[3],1)
    dummy3.missing_x[findall(x->x==entity_idx_no234U238U[3][i],df.entity_id)] =
        dummy4.missing_x[findall(x->x==entity_idx_no234U238U[3][i],df.entity_id)]
end
df.a234U238U = dummy3.missing_x
dummy3 = check_entity(df, entity_idx_no230Th238U, "230Th/238U")
df.a230Th238U = dummy3.missing_x

dummy3 = check_entity(df, entity_idx_no230Th232Thini[1], "a230Th232Thini")
dummy4 = check_entity(df, entity_idx_no230Th232Th[2], "no a230Th/232Th but atomic ratio(?)")
for i = 1:size(entity_idx_no230Th232Thini[2],1)
    dummy3.missing_x[findall(x->x==entity_idx_no230Th232Thini[2][i],df.entity_id)] =
        dummy4.missing_x[findall(x->x==entity_idx_no230Th232Thini[2][i],df.entity_id)]
end
df.a230Th232Thini = dummy3.missing_x

### output of requested csv files (for Laia, sent 30.09.2019),
### wich was later called "datasets_metadata_issues_approval check.csv"
#CSV.write("dataset_details.csv",df)

df2 = DataFrame()
df2 = datingDO[!,[:dating_id,
    :entity_id,:date_type,:depth_dating,:uncorr_age,:uncorr_age_uncert_pos,:uncorr_age_uncert_neg,
    :c238U_content,:c238U_uncertainty,
    :c232Th_content,:c232Th_uncertainty,:c230Th_content,:c230Th_uncertainty,
    :a230Th_232Th_ratio,:a230Th_232Th_ratio_uncertainty,:a230Th_238U_activity,
    :a230Th_238U_activity_uncertainty,:a234U_238U_activity,:a234U_238U_activity_uncertainty,
    :ini_230Th_232Th_ratio,:ini_230Th_232Th_ratio_uncertainty,:decay_constant,
    :corr_age,:corr_age_uncert_pos,:corr_age_uncert_neg]]

### output of requested csv files (for Laia, sent 30.09.2019)
#CSV.write("dating_datasets.csv",df2)

################################################################################
### real quality control of individual dating_ids
### 1. check for all necessary activity ratios or calculate those
################################################################################

k=zeros(Int,length(datingDO.dating_id),7)
for m = 1:length(datingDO.dating_id)
    if ismissing(datingDO.c238U_content[m])
        global k[m,1] = 1
    end
    if ismissing(datingDO.c232Th_content[m])
        global k[m,2] = 1
    end
    if ismissing(datingDO.c230Th_content[m])
        global k[m,3] = 1
    end
    if ismissing(datingDO.a230Th_232Th_ratio[m])
        global k[m,4] = 1
    end
    if ismissing(datingDO.a230Th_238U_activity[m])
        global k[m,5] = 1
    end
    if ismissing(datingDO.a234U_238U_activity[m])
        global k[m,6] = 1
    end
    if ismissing(datingDO.ini_230Th_232Th_ratio[m])
        global k[m,7] = 1
    end
end


################################################################################
# this loop is for preparation of activity ratios (accounting for various data
# availability of ratios and concentrations) and subsequent age calculations
# WITHOUT detrital correction
################################################################################

w1=w2=w4=w6=0
t_uncorr = zeros(length(datingDO.dating_id))
dt_uncorr = zeros(length(datingDO.dating_id))
flag_uncorr = ones(length(datingDO.dating_id))

for m = 1:length(datingDO.dating_id)

    ###############
    # use appropriate decay constants
    ###############

    if datingDO.decay_constant[m] == "Cheng et al. 2000"
        l230 = l230old
        l234 = l234old
    elseif datingDO.decay_constant[m] == "Cheng et al. 2013"
        l230 = l230new
        l234 = l234new
    elseif datingDO.decay_constant[m] == "Edwards et al. 1987"
        l230 = l230veryold
        l234 = l234veryold
    else
        l230 = l230veryold          # necessary to think about that in more detail
        l234 = l234veryold          # necessary to think about that in more detail
    end

    ############################################################################
    # prepare necessary values and calculate age if 0/8 and 4/8 are available
    ############################################################################
    if k[m,5] == 0 && k[m,6] == 0

        UU, dUU = prepare1(datingDO.a234U_238U_activity[m],datingDO.a234U_238U_activity_uncertainty[m],0.01)
        Th0U, dTh0U = prepare1(datingDO.a230Th_238U_activity[m],datingDO.a230Th_238U_activity_uncertainty[m],0.01)

        global w1 += 1

        t=age_err_determination(UU, dUU, Th0U, dTh0U, l230, l234)
        if t[1] != "out of range"
            global t_uncorr[m] = t[1]  # age [a]
            global dt_uncorr[m] = t[2] # 2 sigma error
        end
        #println("1 ",datingDO.dating_id[m]," ", datingDO.uncorr_age[m]," ",t[1]," ",t[2],m)

    ############################################################################
    # prepare necessary values and calculate age if 8, 0 and 4/8 are available
    ############################################################################
    elseif k[m,1] == 0 && k[m,3] == 0 && k[m,6] == 0

        UU, dUU = prepare1(datingDO.a234U_238U_activity[m],datingDO.a234U_238U_activity_uncertainty[m],0.01)
        Th0U, dTh0U = prepare2(datingDO.c230Th_content[m]*l230,datingDO.c230Th_uncertainty[m]*l230,0.01,
            datingDO.c238U_content[m]*l238,datingDO.c238U_uncertainty[m]*l238,0.01)

        global w2 += 1

        t=age_err_determination(UU, dUU, Th0U, dTh0U, l230, l234)
        if t[1] != "out of range"
            global t_uncorr[m] = t[1]  # age [a]
            global dt_uncorr[m] = t[2] # 2 sigma error
        end

    ############################################################################
    # prepare necessary values and calculate age if 8, 2, 0/2 and 4/8 are available
    ############################################################################
    elseif k[m,1] == 0 && k[m,2] == 0 && k[m,4] == 0 && k[m,6] == 0 && k[m,7] == 0

        UU, dUU = prepare1(datingDO.a234U_238U_activity[m],datingDO.a234U_238U_activity_uncertainty[m],0.01)

        Th2U, dTh2U = prepare2(datingDO.c232Th_content[m]*l232,datingDO.c232Th_uncertainty[m]*l232, 0.01,
            datingDO.c238U_content[m]*l238,datingDO.c238U_uncertainty[m]*l238, 0.01)
        Th0U, dTh0U = prepare2(datingDO.a230Th_232Th_ratio[m],datingDO.a230Th_232Th_ratio_uncertainty[m], 0.01,
            1/Th2U, dTh2U/Th2U^2, 0.01/Th2U^2 )

        global w4 += 1

        t=age_err_determination(UU, dUU, Th0U, dTh0U, l230, l234)
        if t[1] != "out of range"
            global t_uncorr[m] = t[1]  # age [a]
            global dt_uncorr[m] = t[2] # 2 sigma error
        end

    else
        #println("no age can be calulated as necessary activity ratios or contents are missing for dating_id $(datingDO.dating_id[m])")
        global w6 += 1
        global flag_uncorr[m] = 0

    end
end

println(w1," ",w2," ",w4," ",w6)
################################################################################


################################################################################
# this loop is for preparation of activity ratios (accounting for various data
# availability of ratios and concentrations) and subsequent age calculations
# INCLUDING detrital correction
################################################################################

m1=m2=m3=m4=m5=m6=0
t_corr = zeros(length(datingDO.dating_id))      #saves corrected ages
dt_corr = zeros(length(datingDO.dating_id))     #saves ages errors of corrected ages
flag_corr = ones(length(datingDO.dating_id))    #if this is ==0; no age could be calculated as some data are missing

for m = 1:length(datingDO.dating_id)

    ###############
    # use appropriate decay constants
    ###############

    if datingDO.decay_constant[m] == "Cheng et al. 2000"
        l230 = l230old
        l234 = l234old
    elseif datingDO.decay_constant[m] == "Cheng et al. 2013"
        l230 = l230new
        l234 = l234new
    else
        l230 = l230old          # necessary to think about that in more detail
        l234 = l234old          # necessary to think about that in more detail
    end

    ############################################################################
    # prepare necessary values and calculate age if 8, 2, 0/8, 4/8 and
    # ini0/2 are available
    ############################################################################
    if k[m,1] == 0 && k[m,2] == 0 && k[m,5] == 0 && k[m,6] == 0 && k[m,7] == 0

        UU, dUU = prepare1(datingDO.a234U_238U_activity[m],datingDO.a234U_238U_activity_uncertainty[m],0.01)
        Th0U, dTh0U = prepare1(datingDO.a230Th_238U_activity[m],datingDO.a230Th_238U_activity_uncertainty[m],0.01)

        Th2U, dTh2U = prepare2(datingDO.c232Th_content[m]*l232/1000,datingDO.c232Th_uncertainty[m]*l232/1000, 0.01,
            datingDO.c238U_content[m]*l238,datingDO.c238U_uncertainty[m]*l238, 0.01)
        Th0Th2ini, dTh0Th2ini = prepare1(datingDO.ini_230Th_232Th_ratio[m],datingDO.ini_230Th_232Th_ratio_uncertainty[m],0.5)

        global m1 += 1

        t=detritus_corr_err(UU, dUU, Th0U, dTh0U, Th2U, dTh2U, Th0Th2ini, dTh0Th2ini, l230, l234)
        if t[1] != "out of range"
            t_corr[m] = t[1]  # age [a]
            dt_corr[m] = t[2] # 2 sigma error
        end
        #println(datingDO.dating_id[m]," ", datingDO.uncorr_age[m]," ",t)

    #######################################################################
    # prepare necessary values and calculate age if 8, 2, 0, 4/8 and
    # ini0/2 are available
    ########################################################################
    elseif k[m,1] == 0 && k[m,2] == 0 && k[m,3] == 0 && k[m,6] == 0 && k[m,7] == 0

        UU, dUU = prepare1(datingDO.a234U_238U_activity[m],datingDO.a234U_238U_activity_uncertainty[m],0.01)
        Th0U, dTh0U = prepare2(datingDO.c230Th_content[m]*l230,datingDO.c230Th_uncertainty[m]*l230,0.01,
            datingDO.c238U_content[m]*l238,datingDO.c238U_uncertainty[m]*l238,0.01)

        Th2U, dTh2U = prepare2(datingDO.c232Th_content[m]*l232/1000,datingDO.c232Th_uncertainty[m]*l232/1000, 0.01,
            datingDO.c238U_content[m]*l238,datingDO.c238U_uncertainty[m]*l238, 0.01)
        Th0Th2ini, dTh0Th2ini = prepare1(datingDO.ini_230Th_232Th_ratio[m],datingDO.ini_230Th_232Th_ratio_uncertainty[m],0.5)

        global m2 += 1

        t=detritus_corr_err(UU, dUU, Th0U, dTh0U, Th2U, dTh2U, Th0Th2ini, dTh0Th2ini, l230, l234)
        if t[1] != "out of range"
            t_corr[m] = t[1]  # age [a]
            dt_corr[m] = t[2] # 2 sigma error
        end
        #println(datingDO.dating_id[m]," ", datingDO.uncorr_age[m]," ",t)

    #######################################################################
    # prepare necessary values and calculate age if 0/2, 0/8, 4/8 and
    # ini0/2 are available
    ########################################################################
    elseif k[m,4] == 0 && k[m,5] == 0 && k[m,6] == 0 && k[m,7] == 0

        UU, dUU = prepare1(datingDO.a234U_238U_activity[m],datingDO.a234U_238U_activity_uncertainty[m],0.01)
        Th0U, dTh0U = prepare1(datingDO.a230Th_238U_activity[m],datingDO.a230Th_238U_activity_uncertainty[m],0.01)

        Th2U, dTh2U = prepare2(datingDO.a230Th_238U_activity[m],datingDO.a230Th_238U_activity_uncertainty[m],0.01,
            datingDO.a230Th_232Th_ratio[m],datingDO.a230Th_232Th_ratio_uncertainty[m], 0.01)
        Th0Th2ini, dTh0Th2ini = prepare1(datingDO.ini_230Th_232Th_ratio[m],datingDO.ini_230Th_232Th_ratio_uncertainty[m],0.5)

        global m3 += 1

        t=detritus_corr_err(UU, dUU, Th0U, dTh0U, Th2U, dTh2U, Th0Th2ini, dTh0Th2ini, l230, l234)
        if t[1] != "out of range"
            t_corr[m] = t[1]  # age [a]
            dt_corr[m] = t[2] # 2 sigma error
        end
        #println(datingDO.dating_id[m]," ", datingDO.uncorr_age[m]," ",t)

    #######################################################################
    # prepare necessary values and calculate age if 8, 2, 0/2, 4/8 and
    # ini0/2 are available
    ########################################################################
    elseif k[m,1] == 0 && k[m,2] == 0 && k[m,4] == 0 && k[m,6] == 0 && k[m,7] == 0

        UU, dUU = prepare1(datingDO.a234U_238U_activity[m],datingDO.a234U_238U_activity_uncertainty[m],0.01)

        Th2U, dTh2U = prepare2(datingDO.c232Th_content[m]*l232/1000,datingDO.c232Th_uncertainty[m]*l232/1000, 0.01,
            datingDO.c238U_content[m]*l238,datingDO.c238U_uncertainty[m]*l238, 0.01)
        Th0U, dTh0U = prepare2(datingDO.a230Th_232Th_ratio[m],datingDO.a230Th_232Th_ratio_uncertainty[m], 0.01,
            1/Th2U, dTh2U/Th2U^2, 0.01/Th2U^2 )

        Th0Th2ini, dTh0Th2ini = prepare1(datingDO.ini_230Th_232Th_ratio[m],datingDO.ini_230Th_232Th_ratio_uncertainty[m],0.5)

        global m4 += 1

        t=detritus_corr_err(UU, dUU, Th0U, dTh0U, Th2U, dTh2U, Th0Th2ini, dTh0Th2ini, l230, l234)
        if t[1] != "out of range"
            t_corr[m] = t[1]  # age [a]
            dt_corr[m] = t[2] # 2 sigma error
        end
        #println(datingDO.dating_id[m]," ", datingDO.uncorr_age[m]," ",t)

    #######################################################################
    # prepare necessary values and calculate age if 2, 0, 0/8, 4/8 and
    # ini0/2 are available
    ########################################################################
    elseif k[m,2] == [0] && k[m,3] == 0 && k[m,5] == 0 && k[m,6] == 0 && k[m,7] == 0

        UU, dUU = prepare1(datingDO.a234U_238U_activity[m],datingDO.a234U_238U_activity_uncertainty[m],0.01)
        Th0U, dTh0U = prepare1(datingDO.a230Th_238U_activity[m],datingDO.a230Th_238U_activity_uncertainty[m],0.01)

        Th2Th0, dTh2Th0 = prepare2(datingDO.c230Th_content[m]*l230,datingDO.c230Th_uncertainty[m]*l230,0.01,
            datingDO.c232Th_content[m]*l232,datingDO.c232Th_uncertainty[m]*l232, 0.01,)
        Th2U, dTh2U = prepare2(datingDO.a230Th_238U_activity[m],datingDO.a230Th_238U_activity_uncertainty[m],0.01,
            1/Th2Th0, dTh2Th0/Th2Th0^2, 0.01/Th2Th0^2)
        Th0Th2ini, dTh0Th2ini = prepare1(datingDO.ini_230Th_232Th_ratio[m],datingDO.ini_230Th_232Th_ratio_uncertainty[m],0.5)

        global m5 += 1

        t=detritus_corr_err(UU, dUU, Th0U, dTh0U, Th2U, dTh2U, Th0Th2ini, dTh0Th2ini, l230, l234)
        if t[1] != "out of range"
            t_corr[m] = t[1]  # age [a]
            dt_corr[m] = t[2] # 2 sigma error
        end
        #println(datingDO.dating_id[m]," ", datingDO.uncorr_age[m]," ",t)

    else
        global m6 += 1
        global flag_corr[m]= 0

    end
end

println(m1," ",m2," ",m3," ",m4," ",m5," ",m6)

########################################################################
# Quality Check of provided data
########################################################################

# find different stalagmites
i = findall(x->x!=0,diff(datingDO.entity_id)).+1
i1= zeros(length(i)+2)
i1[1]=1
i1[2:length(i)+1] = i
i1[length(i)+2] = length(datingDO.entity_id)+1
i1=convert(Vector{Int},i1)      # i1 provides first index of a stal (except of i1(length(i1)) )


# 1. Check for missing values for uncorrected age calculations
for m = 1:length(i1)-1

    ave_flag_uncorr = mean(flag_uncorr[i1[m]: i1[m+1]-1])

    if ave_flag_uncorr == 0
        println("no age can be calulated as necessary activity ratios or contents are missing for dating_ids $(datingDO.dating_id[i1[m]]) - $(datingDO.dating_id[i1[m+1]-1])")
        println("The following table provides an overview, with which set of data it is possible to calculate uncorrected ages.")
        println("number | 238U  |  232Th  |  230Th  |  230Th232Th | 230Th238U | 234U/238U")
        println("  1    |       |         |         |             |     x     |     x    ")
        println("  2    |   x   |         |    x    |             |           |     x    ")
        println("  4    |   x   |    x    |         |      x      |           |     x    ")
        println(" ")
    elseif 0 < ave_flag_uncorr < 1
        println("This is weird. For some of your data datapoints at least one crucial activty ratio is missing")
        println("Please, check for the completeness of your data set. dating_ids: $(datingDO.dating_id[i1[m]]) - $(datingDO.dating_id[i1[m+1]-1])")
        println(" ")
    else
        println("Great, you provided all necessary data to calculate uncorrected ages for dating_ids $(datingDO.dating_id[i1[m]]) - $(datingDO.dating_id[i1[m+1]-1])")

        ave_offset_flag_uncorr1 = abs(mean(t_uncorr[i1[m]:i1[m+1]-1]-datingDO.uncorr_age[i1[m]:i1[m+1]-1]))/
            mean(t_uncorr[i1[m]:i1[m+1]-1])*100 # in %
        println(ave_offset_flag_uncorr1," ",m)
        ave_a234U238U = mean(datingDO.a234U_238U_activity[i1[m]:i1[m+1]-1])

        if ismissing(ave_a234U238U)
            println("not for all ages the activity ratio of 234U/238U is given. Could not check given values.")
            println("")
        elseif ave_a234U238U<0
            println("You have provided d234U instead of activity ratio 234U/238U.")
            println("")
        elseif ave_a234U238U>10
            println("Most likely, you have provided d234U instead of activity ratio 234U/238U.")
            println("")
        elseif ismissing(ave_offset_flag_uncorr1)
            println("not all uncorrected ages provided: not possible to check the given concentrations and activity ratios.")
            println("")
        elseif ave_offset_flag_uncorr1 > 1
            println("Unfortunately, at least one kind of your provided activity ratios or concentrations is not correct.")
            println("Out of experience: the error is most often in activity ratio of 234U/238U (either expressed as d234U or even the initial 234U/238U).")
            println("Alternatively check the other ratios and concentrations and the provided half-lives.")
            println("")
        else
            println("Great, it seems all provided activity ratios to calculate uncorrected ages appear to be correctly included.")
            println("")
        end
    end

    ave_flag_corr = mean(flag_corr[i1[m]: i1[m+1]-1])

    if flag_corr[i1[m]] == 0
        println("no detrial corrected age can be calulated as necessary activity ratios or contents are missing for dating_ids $(datingDO.dating_id[i1[m]]) - $(datingDO.dating_id[i1[m+1]-1])")
        println("The following table provides an overview, with which set of data it is possible to calculate ages.")
        println("number | 238U  |  232Th  |  230Th  |  230Th232Th | 230Th238U | 234U/238U | ini230Th232Th")
        println("  1    |   x   |    x    |         |             |     x     |     x     |       x")
        println("  2    |   x   |    x    |    x    |             |           |     x     |       x")
        println("  3    |       |         |         |      x      |     x     |     x     |       x")
        println("  4    |   x   |    x    |         |      x      |           |     x     |       x")
        println("  5    |       |    x    |    x    |             |     x     |     x     |       x")
        println(" ")
    elseif 0 < ave_flag_corr < 1
        println("This is weird. For some of your data datapoints at least one crucial activty ratio is missing")
        println("Please, check for the completeness of your data set. dating_ids: $(datingDO.dating_id[i1[m]]) - $(datingDO.dating_id[i1[m+1]-1])")
        println(" ")
    else
        println("Great, you provided all necessary data to calculate uncorrected ages for dating_ids $(datingDO.dating_id[i1[m]]) - $(datingDO.dating_id[i1[m+1]-1])")

        ave_offset_flag_corr1 = abs(mean(t_corr[i1[m]:i1[m+1]-1]-datingDO.corr_age[i1[m]:i1[m+1]-1]))/
            mean(t_corr[i1[m]:i1[m+1]-1])*100 # in %
        println(ave_offset_flag_corr1," ",m)
        ave_ini230Th232Th = mean(datingDO.ini_230Th_232Th_ratio[i1[m]:i1[m+1]-1])
        if k[i1[m],3] == 0 || k[i1[m],4] == 0
            if k[i1[m],1] == 0 && k[i1[m],2] == 0
            else
                ave_230Th232Th = mean(datingDO.a230Th_232Th_ratio[i1[m]:i1[m+1]-1])
                println(ave_230Th232Th)
            end
        end

        if ismissing(ave_ini230Th232Th)
            println("not for all ages the initial activity ratio of 230Th/232Th is given. Could not check given values.")
            println("")
        elseif ave_ini230Th232Th<0.01
            println("You have provided initial atomic ratio instead of initial 230Th/232Th activity ratio.")
            println("")
        elseif @isdefined(ave_230Th232Th) && ave_230Th232Th<0.01
            println("You have provided atomic ratio instead of 230Th/232Th activity ratio.")
            println("")
        elseif ismissing(ave_offset_flag_corr1)
            println("not all uncorrected ages provided: not possible to check the given concentrations and activity ratios.")
            println("")
        elseif ave_offset_flag_corr1 > 5
            println("Unfortunately, at least one kind of your provided activity ratios or concentrations is not correct.")
            println("Out of experience: the error is most often in initial activity ratio of 230Th/232Th (either expressed as atomic ratio (+forgotten 10e-6) or as even 232Th/238U atomic ratio).")
            println("Alternatively check the other ratios and concentrations and the provided half-lives.")
            println("")
        else
            println("Great, it seems all provided activity ratios to calculate corrected ages appear to be correctly included.")
            println("")
        end
    end
end
################################################################################

### 2
