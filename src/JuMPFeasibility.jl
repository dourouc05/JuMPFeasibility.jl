module JuMPFeasibility

using LinearAlgebra
using SparseArrays

using JuMP

const _AFF_OR_SV = Union{MOI.ScalarAffineFunction{T}, MOI.SingleVariable} where T
const _LIN_CON = Union{MOI.GreaterThan{T}, MOI.LessThan{T}, MOI.EqualTo{T}} where T

function _lhs_rhs_val_constraint(con_ref::ConstraintRef{Model, MOI.ConstraintIndex{F, S}}, var_value::Function) where {F <: _AFF_OR_SV, S <: _LIN_CON}
	lhs = value(con_ref, var_value)
    rhs = MOI.constant(moi_set(constraint_object(con_ref)))
    return lhs, rhs
end

"""
    is_satisfied(con_ref::ConstraintRef{Model, MOI.ConstraintIndex{F, S}}, var_value::Function; ε::Float64=1.0e-6) where {F, S}

Checks whether the constraint is satisfied by the supplied solution. This simply imples that 
the solution is feasible for this constraint. 
"""
function is_satisfied(con_ref::ConstraintRef{Model, MOI.ConstraintIndex{F, MOI.GreaterThan{T}}}, var_value::Function; ε::Float64=1.0e-6) where {F <: _AFF_OR_SV, T}
	lhs, rhs = _lhs_rhs_val_constraint(con_ref, var_value)
	return lhs >= rhs - ε
end

function is_satisfied(con_ref::ConstraintRef{Model, MOI.ConstraintIndex{F, MOI.LessThan{T}}}, var_value::Function; ε::Float64=1.0e-6) where {F <: _AFF_OR_SV, T}
	lhs, rhs = _lhs_rhs_val_constraint(con_ref, var_value)
	return lhs <= rhs + ε
end

function is_satisfied(con_ref::ConstraintRef{Model, MOI.ConstraintIndex{F, MOI.EqualTo{T}}}, var_value::Function; ε::Float64=1.0e-6) where {F <: _AFF_OR_SV, T}
	lhs, rhs = _lhs_rhs_val_constraint(con_ref, var_value)
	return abs(lhs - rhs) <= ε
end

function is_satisfied(con_ref::ConstraintRef{Model, MOI.ConstraintIndex{F, MOI.Interval{T}}}, var_value::Function; ε::Float64=1.0e-6) where {F <: _AFF_OR_SV, T}
    lhs = value(con_ref, var_value)
    lb = moi_set(constraint_object(con_ref)).lower
    ub = moi_set(constraint_object(con_ref)).upper
	return lb - ε <= lhs <= ub + ε
end

"""
    is_tight(con_ref::ConstraintRef{Model, MOI.ConstraintIndex{F, S}}, var_value::Function; ε::Float64=1.0e-6) where {F, S}

Checks whether the constraint is tight (saturated) at the supplied solution. For an inequality
constraint, this implies that the supplied solution satisfies the constraint at equality. 
"""
function is_tight(con_ref::ConstraintRef{Model, MOI.ConstraintIndex{F, S}}, var_value::Function; ε::Float64=1.0e-6) where {F <: _AFF_OR_SV, S <: _LIN_CON}
	lhs, rhs = _lhs_rhs_val_constraint(con_ref, var_value)
	return abs(lhs - rhs) <= ε 
end

function is_tight(con_ref::ConstraintRef{Model, MOI.ConstraintIndex{F, MOI.Interval{T}}}, var_value::Function; ε::Float64=1.0e-6) where {F <: _AFF_OR_SV, T}
    # An interval constraint cannot be tight, unless it is degenerate (i.e. reduces to one possible value).

    # Check for degeneracy. 
    lb = moi_set(constraint_object(con_ref)).lower
    ub = moi_set(constraint_object(con_ref)).upper
    if abs(ub - lb) > ε
        return false
    end
    
    # As lb ≈ ub, check whether equality is tight.
    lhs = value(con_ref, var_value)
	return min(abs(lhs - ub), abs(lhs - lb)) <= ε
end

"""
    coefficients(con_ref::ConstraintRef{Model, MOI.ConstraintIndex{F, S}}, var_to_idx::Dict{VariableRef, Int}) where {F, S}

Returns the coefficients of the linear constraint in a Julia-native vector-like container. 
The indices are given in argument. If `var_to_idx` lacks some variables, they will not be exported
in the coefficient vector.

TODO: add a function to create this dictionary, with all variables?
"""
function coefficients(con_ref::ConstraintRef{Model, MOI.ConstraintIndex{MOI.ScalarAffineFunction{T}, S}}, var_to_idx::Dict{VariableRef, Int}) where {S, T}
    m = owner_model(con_ref)
    nvars = MOI.get(m, MOI.NumberOfVariables())

    v = spzeros(nvars)
    for (var, coeff) in jump_function(constraint_object(con_ref)).terms
        if var in keys(var_to_idx)
            v[var_to_idx[var]] = coeff
        end
    end
    return v
end

export is_satisfied, is_tight, coefficients

end # module
