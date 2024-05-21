#########################################################################
# shapley_table.jl                                                      #
#                                                                       #
# Makes table of important features based on Shapley indices.           # #                                                                       #
#                                                                       #
# This script requires the Shapley output to be present in              #
#   `output/shapley` .                                                  #
#                                                                       #
#########################################################################

# load environment and packages
import Pkg
Pkg.activate(".")
Pkg.instantiate()

using CSVFiles # read CSV of Shapley indices
using DataFrames # data structure for indices
using PrettyTables # table formatting
using Crayons # crayons for highlighting

# define highlight color for important parameters
crayon_red_light = Crayon(foreground = :white, background = :light_red);

highlight_max = PrettyTables.LatexHighlighter(
  (data, i, j) -> (typeof(data[i, j]) == Float64) && (data[i,j] > 0.10), 
  ["color{red}", "textbf"])

# assume this is called from the project root directory
output_path = joinpath(@__DIR__, "..", "output", "shapley")
shap_ind_def = DataFrame(CSVFiles.load(joinpath(output_path, "shapley_indices_default.csv")))
shap_ind_opt = DataFrame(CSVFiles.load(joinpath(output_path, "shapley_indices_optimistic.csv")))
shap_ind_pes = DataFrame(CSVFiles.load(joinpath(output_path, "shapley_indices_pessimistic.csv")))

# filter indices for most important values
select!(shap_ind_def, [:feature_name, :mean_2050, :mean_2100, :mean_2150, :mean_2200]) # select subseto f columns
# normalize
shap_ind_def[:, 2:end] = reduce(hcat, map(col -> col ./ sum(col), eachcol(shap_ind_def[:, 2:end])))
shap_ind_sig_idx = findall(map(maximum, eachrow(shap_ind_def[:, 2:end])) .> 0.05)
shap_ind_def = shap_ind_def[shap_ind_sig_idx, :]

select!(shap_ind_opt, [:feature_name, :mean_2050, :mean_2100, :mean_2150, :mean_2200]) # select subseto f columns
# normalize
shap_ind_opt[:, 2:end] = reduce(hcat, map(col -> col ./ sum(col), eachcol(shap_ind_opt[:, 2:end])))
shap_ind_sig_idx = findall(map(maximum, eachrow(shap_ind_opt[:, 2:end])) .> 0.05)
shap_ind_opt = shap_ind_opt[shap_ind_sig_idx, :]

select!(shap_ind_pes, [:feature_name, :mean_2050, :mean_2100, :mean_2150, :mean_2200]) # select subseto f columns
# normalize
shap_ind_pes[:, 2:end] = reduce(hcat, map(col -> col ./ sum(col), eachcol(shap_ind_pes[:, 2:end])))
shap_ind_sig_idx = findall(map(maximum, eachrow(shap_ind_pes[:, 2:end])) .> 0.05)
shap_ind_pes = shap_ind_pes[shap_ind_sig_idx, :]

# print LaTeX markup for tables
PrettyTables.pretty_table(shap_ind_def, backend=Val(:latex); highlighters=highlight_max)

PrettyTables.pretty_table(shap_ind_opt, backend=Val(:latex); highlighters=highlight_max)

PrettyTables.pretty_table(shap_ind_pes, backend=Val(:latex); highlighters=highlight_max)
