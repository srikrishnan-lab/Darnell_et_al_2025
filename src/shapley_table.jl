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
output_path = "output/shapley"
shap_ind = DataFrame(CSVFiles.load(joinpath(output_path, "shapley_indices.csv")))
select!(shap_ind, [:feature_name, :mean_2050, :mean_2100, :mean_2150, :mean_2200]) # select subseto f columns
# normalize
shap_ind[:, 2:end] = reduce(hcat, map(col -> col ./ sum(col), eachcol(shap_ind[:, 2:end])))
shap_ind_sig_idx = findall(map(maximum, eachrow(shap_ind[:, 2:end])) .> 0.05)
shap_ind = shap_ind[shap_ind_sig_idx, :]

# print LaTeX markup for table
PrettyTables.pretty_table(shap_ind, backend=Val(:latex); highlighters=highlight_max)