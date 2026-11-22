# this script processes the link between abundances of mammal prey and body size
# this is broken up into a Bayesian section and a GLM section, which will use different packages

##### specify folder and load data
dir = "C:/Users/93652672/OneDrive - University of Nebraska/Projects in progress/Barn_owl_global_diet/"
cd(dir)

# need these for all parts
using CSV
using DataFrames
using Statistics
using CategoricalArrays
using GLMakie
using GLM

using CairoMakie

colorgrad = cgrad(:roma, 15, categorical = true, scale = :exp)

# read in the data
df_density = CSV.read("Dataset_density2.csv",DataFrame) # read in density data
df_density = df_density[findall(df_density.Order .!= "NA"),:] # drop everything that is rows with just NA
#df_dens_rodents = df_density[findall(df_density.Order .== "Rodentia" .|| df_density.Order .== "Eulipotyphla" .|| df_density.Order .== "Lagomorpha" .|| df_density.Order .== "Dasyuromorphia" .|| df_density.Order .== "Didelphimorphia"),:] # drop everything that is rows with just NA

number_rows = size(df_density,1)
masses = zeros(Float64,number_rows)

df_bodymass = CSV.read("Dataset_predictions.csv",DataFrame) # read in density data

# litte loop to pull mammal body masses from the Elton database and match them to the density database
for i = 1:number_rows
    row_species = df_density.AcceptedName_COL[i]
    mass_to_fill = df_bodymass.BM[findall(df_bodymass.Species .== row_species)]
    if length(mass_to_fill) > 1
        mass_to_fill = mean(mass_to_fill)
    end
    if !isempty(mass_to_fill)
        masses[i] = mass_to_fill[]
    end
end

# add the body mass and logged density data to the main database
df_density[!, :BodyMass] .= masses
df_density[!, :LogMass] = log10.(df_density.BodyMass) # 
df_density[!, :LogDensity] = log10.(df_density.Density_km) # 
df_density = df_density[findall(df_density.BodyMass .!= 0),:] # drop rows with 0 body mass


# write a little function to make sure that the predicitons are floats
function setpredictiontypes(df_predict_y) 
    df_predict_y[!,:prediction] = convert.(Float64,df_predict_y[!,:prediction])
    df_predict_y[!,:lower] = convert.(Float64,df_predict_y[!,:lower])
    df_predict_y[!,:upper] = convert.(Float64,df_predict_y[!,:upper])
end

################################################################
######## set up to do density allometry with a OLS approach
################################################################

########### first do standard Damuth's rule ###############################

f_density = Figure()
ax_all = Axis(f_density[1, 1], xlabel = "log10 Body mass (g)", ylabel = "log10 Density (individuals/km2)") #,yscale = log10,xscale = log10
scatter!(ax_all,df_density.LogMass,df_density.LogDensity, color = colorgrad[8], alpha = 0.6)

lm_all = lm(@formula(LogDensity ~ LogMass), df_density)
allr2 = r2(lm_all)

# make model predictions and plot
new_x = DataFrame(LogMass = range(minimum(df_density.LogMass), maximum(df_density.LogMass), length=10))
df_predict_y = predict(lm_all,new_x, interval=:confidence, level = 0.95)
# ensure that all data types are floats with our above function
setpredictiontypes(df_predict_y)
band!(ax_all,new_x.LogMass,df_predict_y.lower,df_predict_y.upper; color = colorgrad[8])
rod = lines!(ax_all,new_x.LogMass, df_predict_y.prediction, linewidth = 2, color = :black)

text!(ax_all, 0.2, 4, text = "A")
############# now do it by order ########################

# add a panel to the figure
ax_order = Axis(f_density[2, 1], xlabel = "log10 Body mass (g)", ylabel = "log10 Density (individuals/km2)") #,yscale = log10,xscale = log10

# subset the data by the orders actually in the barn owl diet
subset_df = df_density

# tally up rows by order to check on sample sizes
sums_by_order = combine(groupby(subset_df, :Order), nrow)

prey_orders = ["Eulipotyphla","Rodentia","Didelphimorphia","Chiroptera","Lagomorpha","Carnivora","Dasyuromorphia","Macroscelidea","Diprotodontia","Afrosoricida","Primates","Peramelemorphia","Microbiotheria"]
missing_words = setdiff(subset_df.Order, prey_orders)
for j = 1:length(missing_words)
    subset_df = subset_df[findall(subset_df.Order .!= missing_words[j]),:] # drop rows with orders that don't occur in the diets
end

    subset_df = subset_df[findall(subset_df.Order .!= "Microbiotheria"),:] # drop rows, only two diet items and no densities
    subset_df = subset_df[findall(subset_df.Order .!= "Peramelemorphia"),:] # drop rows with 0 body mass

unique(subset_df.Order)

# fit the linear models
lm_1 = lm(@formula(LogDensity ~ LogMass + Order), subset_df)
mainr2 = r2(lm_1)

lm_2 = lm(@formula(LogDensity ~ LogMass * Order), subset_df)
mainr2 = r2(lm_2)

dof_residual(lm_2)
dof(lm_2)


# plot all the points first, then come back to the lines
# rodents
    df_subset = subset_df[findall(subset_df.Order .== "Rodentia"),:]
    scatter!(ax_order,df_subset.LogMass,df_subset.LogDensity, color = ("White", 0.1), strokewidth = 0.5, strokecolor = colorgrad[1])
# rabbits
    df_subset = subset_df[findall(subset_df.Order .== "Lagomorpha"),:]
    scatter!(ax_order,df_subset.LogMass,df_subset.LogDensity, color = ("White", 0.1), strokewidth = 0.5, strokecolor = colorgrad[3])
# insectivores
    df_subset = subset_df[findall(subset_df.Order .== "Eulipotyphla"),:]
    scatter!(ax_order,df_subset.LogMass,df_subset.LogDensity, color = ("White", 0.1), strokewidth = 0.5, strokecolor = colorgrad[5])
# carnivorous marsupials
    df_subset = subset_df[findall(subset_df.Order .== "Dasyuromorphia"),:]
    scatter!(ax_order,df_subset.LogMass,df_subset.LogDensity, color = ("White", 0.1), strokewidth = 0.5, strokecolor = colorgrad[7])
# carnivora
    df_subset = subset_df[findall(subset_df.Order .== "Carnivora"),:]
    scatter!(ax_order,df_subset.LogMass,df_subset.LogDensity, color = ("White", 0.1), strokewidth = 0.5, strokecolor = colorgrad[9])
# primates
    df_subset = subset_df[findall(subset_df.Order .== "Primates"),:]
    scatter!(ax_order,df_subset.LogMass,df_subset.LogDensity, color = ("White", 0.1), strokewidth = 0.5, strokecolor = colorgrad[11])
# Diprotodontia
    df_subset = subset_df[findall(subset_df.Order .== "Diprotodontia"),:]
    scatter!(ax_order,df_subset.LogMass,df_subset.LogDensity, color = ("White", 0.1), strokewidth = 0.5, strokecolor = colorgrad[13])
# Didelphimorphia
    df_subset = subset_df[findall(subset_df.Order .== "Didelphimorphia"),:]
    scatter!(ax_order,df_subset.LogMass,df_subset.LogDensity, color = ("White", 0.1), strokewidth = 0.5, strokecolor = colorgrad[15])

# now get and plot the fitted lines
    # rodents
        df_subset = subset_df[findall(subset_df.Order .== "Rodentia"),:]
        new_x = DataFrame(LogMass = range(minimum(df_subset.LogMass), maximum(df_subset.LogMass), length=10))
        new_x[!, :Order] = vcat(fill.("Rodentia",10))
        df_predict_y = predict(lm_2,new_x, interval=:confidence, level = 0.95)
        # ensure that all data types are floats with our above function
        setpredictiontypes(df_predict_y)
        band!(ax_order,new_x.LogMass,df_predict_y.lower,df_predict_y.upper; color = colorgrad[1])
        rod = lines!(ax_order,new_x.LogMass, df_predict_y.prediction, linewidth = 2, color = :black)
    # rabbits
        df_subset = subset_df[findall(subset_df.Order .== "Lagomorpha"),:]
        new_x = DataFrame(LogMass = range(minimum(df_subset.LogMass), maximum(df_subset.LogMass), length=10))
        new_x[!, :Order] = vcat(fill.("Lagomorpha",10))
        df_predict_y = predict(lm_2,new_x, interval=:confidence, level = 0.95)
        # ensure that all data types are floats with our above function
        setpredictiontypes(df_predict_y)
        band!(ax_order,new_x.LogMass,df_predict_y.lower,df_predict_y.upper; color = colorgrad[3])
        lag = lines!(ax_order,new_x.LogMass, df_predict_y.prediction, linewidth = 2, color = :black)
    # plot the line for insectivores
        df_subset = subset_df[findall(subset_df.Order .== "Eulipotyphla"),:]
        new_x = DataFrame(LogMass = range(minimum(df_subset.LogMass), maximum(df_subset.LogMass), length=10))
        new_x[!, :Order] = vcat(fill.("Eulipotyphla",10))
        df_predict_y = predict(lm_2,new_x, interval=:confidence, level = 0.95)
        # ensure that all data types are floats with our above function
        setpredictiontypes(df_predict_y)
        band!(ax_order,new_x.LogMass,df_predict_y.lower,df_predict_y.upper; color = colorgrad[5])
        ins = lines!(ax_order,new_x.LogMass, df_predict_y.prediction, linewidth = 2, color = :black)
    # plot the line for carnivorous marsupials
        df_subset = subset_df[findall(subset_df.Order .== "Dasyuromorphia"),:]
        new_x = DataFrame(LogMass = range(minimum(df_subset.LogMass), maximum(df_subset.LogMass), length=10))
        new_x[!, :Order] = vcat(fill.("Dasyuromorphia",10))
        df_predict_y = predict(lm_2,new_x, interval=:confidence, level = 0.95)
        # ensure that all data types are floats with our above function
        setpredictiontypes(df_predict_y)
        band!(ax_order,new_x.LogMass,df_predict_y.lower,df_predict_y.upper; color = colorgrad[7])
        cm = lines!(ax_order,new_x.LogMass, df_predict_y.prediction, linewidth = 2, color = :black)
    # plot the line for carnivora
        df_subset = subset_df[findall(subset_df.Order .== "Carnivora"),:]
        new_x = DataFrame(LogMass = range(minimum(df_subset.LogMass), maximum(df_subset.LogMass), length=10))
        new_x[!, :Order] = vcat(fill.("Carnivora",10))
        df_predict_y = predict(lm_2,new_x, interval=:confidence, level = 0.95)
        # ensure that all data types are floats with our above function
        setpredictiontypes(df_predict_y)
        band!(ax_order,new_x.LogMass,df_predict_y.lower,df_predict_y.upper; color = colorgrad[9])
        cv = lines!(ax_order,new_x.LogMass, df_predict_y.prediction, linewidth = 2, color = :black)
    # plot the line for Primates
        df_subset = subset_df[findall(subset_df.Order .== "Primates"),:]
        new_x = DataFrame(LogMass = range(minimum(df_subset.LogMass), maximum(df_subset.LogMass), length=10))
        new_x[!, :Order] = vcat(fill.("Primates",10))
        df_predict_y = predict(lm_2,new_x, interval=:confidence, level = 0.95)
        # ensure that all data types are floats with our above function
        setpredictiontypes(df_predict_y)
        band!(ax_order,new_x.LogMass,df_predict_y.lower,df_predict_y.upper; color = colorgrad[11])
        pri = lines!(ax_order,new_x.LogMass, df_predict_y.prediction, linewidth = 2, color = :black)
    # plot the line for Diprotodontia
        df_subset = subset_df[findall(subset_df.Order .== "Diprotodontia"),:]
        new_x = DataFrame(LogMass = range(minimum(df_subset.LogMass), maximum(df_subset.LogMass), length=10))
        new_x[!, :Order] = vcat(fill.("Diprotodontia",10))
        df_predict_y = predict(lm_2,new_x, interval=:confidence, level = 0.95)
        # ensure that all data types are floats with our above function
        setpredictiontypes(df_predict_y)
        band!(ax_order,new_x.LogMass,df_predict_y.lower,df_predict_y.upper; color = colorgrad[13])
        dip = lines!(ax_order,new_x.LogMass, df_predict_y.prediction, linewidth = 2, color = :black)
    # plot the line for Didelphimorphia
        df_subset = subset_df[findall(subset_df.Order .== "Didelphimorphia"),:]
        new_x = DataFrame(LogMass = range(minimum(df_subset.LogMass), maximum(df_subset.LogMass), length=10))
        new_x[!, :Order] = vcat(fill.("Didelphimorphia",10))
        df_predict_y = predict(lm_2,new_x, interval=:confidence, level = 0.95)
        # ensure that all data types are floats with our above function
        setpredictiontypes(df_predict_y)
        band!(ax_order,new_x.LogMass,df_predict_y.lower,df_predict_y.upper; color = colorgrad[15])
        did = lines!(ax_order,new_x.LogMass, df_predict_y.prediction, linewidth = 2, color = :black)

text!(ax_order, 0.2, 4, text = "B")

f_density

Legend(f_density[2, 2],
    [PolyElement(color = colorgrad[1]),PolyElement(color = colorgrad[3]),PolyElement(color = colorgrad[5]),
    PolyElement(color = colorgrad[7]),PolyElement(color = colorgrad[9]),PolyElement(color = colorgrad[1]),
    PolyElement(color = colorgrad[13]),PolyElement(color = colorgrad[15])],
    ["Rodentia", "Lagomorphs", "Insectivores", "Carnivorous marsupials", "Carnivora", 
    "Primates", "Diprotodontia", "Didelphimorphia"])

    save("Mammal_density_allometry.pdf",f_density)
