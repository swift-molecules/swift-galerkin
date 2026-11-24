extension Galerkin {
    /// A reference cell's shape functions at one point: each one's value and its derivatives along ξ and η, and
    /// along ζ on a solid cell, one per corner in the cell's own order.
    public struct Shape<Scalar: BinaryFloatingPoint> {
        /// Nᵢ at the point.
        public let values: [Scalar]
        /// ∂Nᵢ/∂ξ at the point.
        public let xi: [Scalar]
        /// ∂Nᵢ/∂η at the point.
        public let eta: [Scalar]
        /// ∂Nᵢ/∂ζ at the point on a solid cell; empty on a flat one.
        public let zeta: [Scalar]

        @inlinable
        public init(values: [Scalar], xi: [Scalar], eta: [Scalar], zeta: [Scalar] = []) {
            precondition(values.count == xi.count && values.count == eta.count, "one value and gradient per corner")
            precondition(zeta.isEmpty || zeta.count == values.count, "no ζ, or one per corner")
            self.values = values
            self.xi = xi
            self.eta = eta
            self.zeta = zeta
        }
    }
}

extension Galerkin.Shape {
    /// How many shape functions there are: one per corner.
    @inlinable
    public var count: Int { values.count }

    /// Σ Nᵢ cᵢ: the interpolation of the corner values `corners`.
    @inlinable
    public func interpolated(_ corners: [Scalar]) -> Scalar {
        precondition(corners.count == count, "one value per corner")
        var sum: Scalar = 0
        for index in values.indices { sum += values[index] * corners[index] }
        return sum
    }
}

extension Galerkin.Shape: Equatable {}
extension Galerkin.Shape: Hashable where Scalar: Hashable {}
extension Galerkin.Shape: Sendable where Scalar: Sendable {}
