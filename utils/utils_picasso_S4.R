setClass("picasso",
    slots = list(
        variables = "character",
        packages = "character",
        work_dir = "character",
        input_path = "character",
        output_path = "character"
    )
)

myObj <- new("picasso",
    variables = ls(), # as.list(environment()),
    packages = .packages(),
    work_dir = "work",
    input_path = "input.txt",
    output_path = "output.txt"
)
