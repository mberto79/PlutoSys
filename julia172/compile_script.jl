using PackageCompiler
# using CSV
# using DataFrames
# using LaTeXStrings
# using Latexify
# using Unitful, UnitfulLatexify, UnitfulRecipes
# using Plots
# using Pluto
# using PlutoUI
# using Markdown

create_sysimage([
    :CSV, 
    :DataFrames,
    :LaTeXStrings,
    :Latexify, :Unitful, :UnitfulLatexify, :UnitfulRecipes,
    :Plots,
    :Pluto,
    :PlutoUI,
    :Markdown
    ];
    precompile_statements_file="pluto_precompile.jl",
    precompile_execution_file="pluto_warmup.jl",
    sysimage_path="PlutoSysImage.so")