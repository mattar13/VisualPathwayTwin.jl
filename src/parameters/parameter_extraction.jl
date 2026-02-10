# ============================================================
# parameter_extraction.jl - Extract parameters from CSV files
# ============================================================

using CSV
using DataFrames

"""
    load_photoreceptor_params_from_csv(csv_path::String)

Load photoreceptor parameters from a CSV file and return as a NamedTuple.
The CSV should have columns: Key, Value, LowerBounds, UpperBounds, DEFAULT

# Arguments
- `csv_path`: Path to the CSV file

# Returns
- Named tuple with parameter names as symbols and values as numbers
"""
function load_photoreceptor_params_from_csv(csv_path::String)
    # Read CSV file
    df = CSV.read(csv_path, DataFrame)

    # Extract parameter names (as symbols) and values
    param_names = Symbol.(df.Key)
    param_values = df.Value

    # Create NamedTuple
    return NamedTuple{Tuple(param_names)}(param_values)
end

"""
    default_rod_params_csv()

Load default rod photoreceptor parameters from the bundled CSV file.
"""
function default_rod_params_csv()
    csv_path = joinpath(@__DIR__, "photoreceptor_params.csv")
    return load_photoreceptor_params_from_csv(csv_path)
end

"""
    default_hc_params_csv()

Load default horizontal cell parameters from the bundled CSV file.
"""
function default_hc_params_csv()
    csv_path = joinpath(@__DIR__, "horizontal_params.csv")
    return load_photoreceptor_params_from_csv(csv_path)
end

"""
    default_on_bc_params_csv()

Load default ON bipolar cell parameters from the bundled CSV file.
"""
function default_on_bc_params_csv()
    csv_path = joinpath(@__DIR__, "on_bipolar_params.csv")
    return load_photoreceptor_params_from_csv(csv_path)
end

"""
    default_off_bc_params_csv()

Load default OFF bipolar cell parameters from the bundled CSV file.
"""
function default_off_bc_params_csv()
    csv_path = joinpath(@__DIR__, "off_bipolar_params.csv")
    return load_photoreceptor_params_from_csv(csv_path)
end

"""
    default_a2_params_csv()

Load default A2 amacrine cell parameters from the bundled CSV file.
"""
function default_a2_params_csv()
    csv_path = joinpath(@__DIR__, "a2_amacrine_params.csv")
    return load_photoreceptor_params_from_csv(csv_path)
end

"""
    default_gaba_params_csv()

Load default GABAergic amacrine cell parameters from the bundled CSV file.
"""
function default_gaba_params_csv()
    csv_path = joinpath(@__DIR__, "gaba_amacrine_params.csv")
    return load_photoreceptor_params_from_csv(csv_path)
end

"""
    default_da_params_csv()

Load default dopaminergic amacrine cell parameters from the bundled CSV file.
"""
function default_da_params_csv()
    csv_path = joinpath(@__DIR__, "da_amacrine_params.csv")
    return load_photoreceptor_params_from_csv(csv_path)
end

"""
    default_gc_params_csv()

Load default ganglion cell parameters from the bundled CSV file.
"""
function default_gc_params_csv()
    csv_path = joinpath(@__DIR__, "ganglion_params.csv")
    return load_photoreceptor_params_from_csv(csv_path)
end

"""
    default_muller_params_csv()

Load default Müller glial cell parameters from the bundled CSV file.
"""
function default_muller_params_csv()
    csv_path = joinpath(@__DIR__, "muller_params.csv")
    return load_photoreceptor_params_from_csv(csv_path)
end

"""
    default_rpe_params_csv()

Load default RPE cell parameters from the bundled CSV file.
"""
function default_rpe_params_csv()
    csv_path = joinpath(@__DIR__, "rpe_params.csv")
    return load_photoreceptor_params_from_csv(csv_path)
end

# ============================================================
# Global parameter and state management
# ============================================================

"""
    load_all_params()

Load all cell type parameters from CSV files and return as a single NamedTuple.

# Returns
A NamedTuple with fields:
- `rod`: Rod photoreceptor parameters
- `hc`: Horizontal cell parameters
- `on_bc`: ON bipolar cell parameters
- `off_bc`: OFF bipolar cell parameters
- `a2`: A2 amacrine cell parameters
- `gaba`: GABAergic amacrine cell parameters
- `da`: Dopaminergic amacrine cell parameters
- `gc`: Ganglion cell parameters
- `muller`: Müller glial cell parameters
- `rpe`: RPE cell parameters

# Example
```julia
params = load_all_params()
rod_C_m = params.rod.C_m
on_bc_g_Ca = params.on_bc.g_Ca
```
"""
function load_all_params()
    return (
        rod = default_rod_params_csv(),
        hc = default_hc_params_csv(),
        on_bc = default_on_bc_params_csv(),
        off_bc = default_off_bc_params_csv(),
        a2 = default_a2_params_csv(),
        gaba = default_gaba_params_csv(),
        da = default_da_params_csv(),
        gc = default_gc_params_csv(),
        muller = default_muller_params_csv(),
        rpe = default_rpe_params_csv()
    )
end
