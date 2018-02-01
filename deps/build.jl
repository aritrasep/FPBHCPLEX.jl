if !("CPLEXExtensions" in keys(Pkg.installed()))
	Pkg.clone("https://github.com/aritrasep/CPLEXExtensions.jl")
	Pkg.build("CPLEXExtensions")
end
