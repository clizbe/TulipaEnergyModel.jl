export resolution_matrix, compute_rp_periods

"""
    M = resolution_matrix(rp_periods, time_steps)

Computes the resolution balance matrix using the array of `rp_periods` and the array of `time_steps`.
The elements in these arrays must be ranges.

## Examples

The following two examples are for two flows/assets with resolutions of 3h and 4h, so that the representative period has 4h periods.

```jldoctest
rp_periods = [1:4, 5:8, 9:12]
time_steps = [1:4, 5:8, 9:12]
resolution_matrix(rp_periods, time_steps)

# output

3×3 Matrix{Float64}:
 1.0  0.0  0.0
 0.0  1.0  0.0
 0.0  0.0  1.0
```

```jldoctest
rp_periods = [1:4, 5:8, 9:12]
time_steps = [1:3, 4:6, 7:9, 10:12]
resolution_matrix(rp_periods, time_steps)

# output

3×4 Matrix{Float64}:
 1.0  0.333333  0.0       0.0
 0.0  0.666667  0.666667  0.0
 0.0  0.0       0.333333  1.0
```
"""
function resolution_matrix(
    rp_periods::AbstractVector{<:UnitRange{<:Integer}},
    time_steps::AbstractVector{<:UnitRange{<:Integer}},
)
    matrix = [
        length(period ∩ time_step) / length(time_step) for period in rp_periods,
        time_step in time_steps
    ]

    return matrix
end

"""
    rp_periods = compute_rp_periods(array_time_steps)

Given the time steps of various flows/assets in the `array_time_steps` input, compute the representative period splits.
Each element of `array_time_steps` is an array of ranges with the following assumptions:

  - An element is of the form `V = [r₁, r₂, …, rₘ]`, where each `rᵢ` is a range `a:b`.
  - `r₁` starts at 1.
  - `rᵢ₊₁` starts at the end of `rᵢ` plus 1.
  - `rₘ` ends at some value `N`, that is the same for all elements of `array_time_steps`.

Notice that this implies that they form a disjunct partition of `1:N`.

The output will also be an array of ranges with the conditions above.
The output is constructed greedily, i.e., it selects the next largest breakpoint following the algorithm below:

 0. Input: `Vᴵ₁, …, Vᴵₚ`, a list of time step ranges. Each element of `Vᴵⱼ` is a range `r = r.start:r.end`. Output: `V`.
 1. Compute the end of the representative period `N` (all `Vᴵⱼ` should have the same end)
 2. Start with an empty `V = []`
 3. Define the beginning of the range `s = 1`
 4. Define an array with all the next breakpoints `B` such that `Bⱼ` is the first `r.end` such that `r.end ≥ s` for each `r ∈ Vᴵⱼ`.
 5. The end of the range will be the `e = max Bⱼ`.
 6. Define `r = s:e` and add `r` to the end of `V`.
 7. If `e = N`, then END
 8. Otherwise, define `s = e + 1` and go to step 4.

## Examples

```jldoctest
time_steps1 = [1:4, 5:8, 9:12]
time_steps2 = [1:3, 4:6, 7:9, 10:12]
compute_rp_periods([time_steps1, time_steps2])

# output

3-element Vector{UnitRange{Int64}}:
 1:4
 5:8
 9:12
```

```jldoctest
time_steps1 = [1:1, 2:3, 4:6, 7:10, 11:12]
time_steps2 = [1:2, 3:4, 5:5, 6:7, 8:9, 10:12]
compute_rp_periods([time_steps1, time_steps2])

# output

5-element Vector{UnitRange{Int64}}:
 1:2
 3:4
 5:6
 7:10
 11:12
```
"""
function compute_rp_periods(
    array_time_steps::AbstractVector{<:AbstractVector{<:UnitRange{<:Integer}}},
)
    # Get Vᴵ₁, the last range of it, the last element of the range
    representative_period_end = array_time_steps[1][end][end]
    for time_steps in array_time_steps
        # Assumption: All start at 1 and end at N
        @assert time_steps[1][1] == 1
        @assert representative_period_end == time_steps[end][end]
    end
    rp_periods = UnitRange{Int}[] # List of ranges
    period_start = 1

    while period_start < representative_period_end
        # The first range end larger than period_start for each range in each time_steps.
        breakpoints = (
            first(r[end] for r in time_steps if r[end] ≥ period_start) for
            time_steps in array_time_steps
        )
        period_end = maximum(breakpoints)
        @assert period_end ≥ period_start
        push!(rp_periods, period_start:period_end)
        period_start = period_end + 1
    end
    return rp_periods
end
