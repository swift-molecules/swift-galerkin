extension Galerkin.Jacobian {
    /// The derivative of the map from a solid reference cell onto an element at one point: the derivatives of
    /// (x, y, z) along ξ, along η and along ζ.
    public struct Solid {
        /// ∂(x, y, z)/∂ξ.
        public let xi: Vector<3, Scalar>
        /// ∂(x, y, z)/∂η.
        public let eta: Vector<3, Scalar>
        /// ∂(x, y, z)/∂ζ.
        public let zeta: Vector<3, Scalar>

        @inlinable
        public init(xi: Vector<3, Scalar>, eta: Vector<3, Scalar>, zeta: Vector<3, Scalar>) {
            self.xi = xi
            self.eta = eta
            self.zeta = zeta
        }
    }
}

extension Galerkin.Jacobian.Solid {
    /// The Jacobian of the element with these corners, in the reference cell's order, where `shape` was evaluated.
    @inlinable
    public init(_ shape: Galerkin.Shape<Scalar>, corners: [Point<3, Scalar>]) {
        precondition(corners.count == shape.count, "one corner per shape function")
        precondition(shape.zeta.count == shape.count, "a solid cell's shape")
        var xi = Vector<3, Scalar>.zero
        var eta = Vector<3, Scalar>.zero
        var zeta = Vector<3, Scalar>.zero
        for index in corners.indices {
            let corner = corners[index].coordinates
            xi += corner.scaled(by: shape.xi[index])
            eta += corner.scaled(by: shape.eta[index])
            zeta += corner.scaled(by: shape.zeta[index])
        }
        self.init(xi: xi, eta: eta, zeta: zeta)
    }

    /// The triple product x,ξ · (x,η × x,ζ): the ratio of an element's volume to the reference cell's at this
    /// point, positive where the map keeps the reference cell's orientation.
    @inlinable
    public var determinant: Scalar { xi.dot(eta.cross(zeta)) }

    /// The derivatives along x, y and z of a quantity whose derivatives along ξ, η and ζ are `xi`, `eta` and
    /// `zeta`: J⁻¹ applied to (∂/∂ξ, ∂/∂η, ∂/∂ζ), whose columns are x,η × x,ζ, x,ζ × x,ξ and x,ξ × x,η over the
    /// determinant.
    @inlinable
    public func cartesian(xi a: Scalar, eta b: Scalar, zeta c: Scalar) -> (x: Scalar, y: Scalar, z: Scalar) {
        let sum = eta.cross(zeta).scaled(by: a) + zeta.cross(xi).scaled(by: b) + xi.cross(eta).scaled(by: c)
        let det = determinant
        return (sum.x / det, sum.y / det, sum.z / det)
    }

    /// Every shape function's derivatives along x, y and z.
    @inlinable
    public func gradients(_ shape: Galerkin.Shape<Scalar>) -> (x: [Scalar], y: [Scalar], z: [Scalar]) {
        precondition(shape.zeta.count == shape.count, "a solid cell's shape")
        var x: [Scalar] = []
        var y: [Scalar] = []
        var z: [Scalar] = []
        x.reserveCapacity(shape.count)
        y.reserveCapacity(shape.count)
        z.reserveCapacity(shape.count)
        for index in 0..<shape.count {
            let gradient = cartesian(xi: shape.xi[index], eta: shape.eta[index], zeta: shape.zeta[index])
            x.append(gradient.x)
            y.append(gradient.y)
            z.append(gradient.z)
        }
        return (x, y, z)
    }
}

extension Galerkin.Jacobian.Solid: Equatable {}
extension Galerkin.Jacobian.Solid: Hashable where Scalar: Hashable {}
extension Galerkin.Jacobian.Solid: Sendable where Scalar: Sendable {}
