
"""
    AbstractBoundaryCondition
格子の境界条件の抽象型。
- PBC: 周期境界条件
- OBC: 開境界条件
- SSD: Sine-Square Deformation 境界条件
"""
abstract type AbstractBoundaryCondition end
struct PBC <: AbstractBoundaryCondition end
struct OBC <: AbstractBoundaryCondition end
export PBC, OBC