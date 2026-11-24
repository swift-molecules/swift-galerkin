extension Galerkin {
    /// The derivative of the map from a reference cell onto an element, [[x,ξ  y,ξ], [x,η  y,η]], at one point.
    public struct Jacobian<Scalar: BinaryFloatingPoint> {
        /// ∂x/∂ξ.
        public let xXi: Scalar
        /// ∂y/∂ξ.
        public let yXi: Scalar
        /// ∂x/∂η.
        public let xEta: Scalar
        /// ∂y/∂η.
        public let yEta: Scalar

        @inlinable
        public init(xXi: Scalar, yXi: Scalar, xEta: Scalar, yEta: Scalar) {
            self.xXi = xXi
            self.yXi = yXi
            self.xEta = xEta
            self.yEta = yEta
        }
    }
}

extension Galerkin.Jacobian {
    /// The Jacobian of the element with these corners, in the reference cell's order, where `shape` was evaluated.
    @inlinable
    public init(_ shape: Galerkin.Shape<Scalar>, corners: [Point<2, Scalar>]) {
        precondition(corners.count == shape.count, "one corner per shape function")
        var xXi: Scalar = 0
        var yXi: Scalar = 0
        var xEta: Scalar = 0
        var yEta: Scalar = 0
        for index in corners.indices {
            xXi += shape.xi[index] * corners[index].x
            yXi += shape.xi[index] * corners[index].y
            xEta += shape.eta[index] * corners[index].x
            yEta += shape.eta[index] * corners[index].y
        }
        self.init(xXi: xXi, yXi: yXi, xEta: xEta, yEta: yEta)
    }

    /// x,ξ y,η − y,ξ x,η: the ratio of an element's area to the reference cell's at this point, positive where the
    /// map keeps the reference cell's orientation.
    @inlinable
    public var determinant: Scalar { xXi * yEta - yXi * xEta }

    /// The derivatives along x and y of a quantity whose derivatives along ξ and η are `xi` and `eta`: J⁻¹ applied
    /// to (∂/∂ξ, ∂/∂η). A covariant vector's components turn into Cartesian ones the same way.
    @inlinable
    public func cartesian(xi: Scalar, eta: Scalar) -> (x: Scalar, y: Scalar) {
        let det = determinant
        return ((yEta * xi - yXi * eta) / det, (xXi * eta - xEta * xi) / det)
    }

    /// Every shape function's derivatives along x and y.
    @inlinable
    public func gradients(_ shape: Galerkin.Shape<Scalar>) -> (x: [Scalar], y: [Scalar]) {
        var x: [Scalar] = []
        var y: [Scalar] = []
        x.reserveCapacity(shape.count)
        y.reserveCapacity(shape.count)
        for index in 0..<shape.count {
            let gradient = cartesian(xi: shape.xi[index], eta: shape.eta[index])
            x.append(gradient.x)
            y.append(gradient.y)
        }
        return (x, y)
    }
}

extension Galerkin.Jacobian: Equatable {}
extension Galerkin.Jacobian: Hashable where Scalar: Hashable {}
extension Galerkin.Jacobian: Sendable where Scalar: Sendable {}
