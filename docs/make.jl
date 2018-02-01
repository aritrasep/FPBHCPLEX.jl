using Documenter, FPBHCPLEX

makedocs(modules=[FPBHCPLEX],
         doctest = false,
         format = :html,
         sitename = "FPBHCPLEX",
         authors = "Aritra Pal",
         pages = Any[
        	"Home" => "index.md",
        	"Installation" => "installation.md",
        	"Getting Started" => "getting_started.md",
        	"Advanced Features" => "advanced.md",
        	"Solving Instances from Literature" => "solving_instances_from_literature.md"
    	])

deploydocs(
	repo = "github.com/aritrasep/FPBHCPLEX.jl.git",
    target = "build",
    osname = "linux",
    julia  = "0.6",
    deps   = nothing,
    make   = nothing,
)
