import pyprop


# Load
input_file = "dev/sample_search.json"
opter = pyprop.Optimizer()
opter.random_search(input_file)