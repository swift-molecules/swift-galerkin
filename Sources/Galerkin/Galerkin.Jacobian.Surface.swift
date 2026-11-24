extension Galerkin.Jacobian {
    /// The derivative of the map from a flat reference cell onto a face in space at one point: the two tangents
    /// ∂(x, y, z)/∂ξ and ∂(x, y, z)/∂η, such as a shell's midsurface or the wall a body of air is bounded by.
    public struct Surface {
        /// ∂(x, y, z)/∂ξ.
        public let xi: Vector<3, Scalar>
        /// ∂(x, y, z)/∂η.
        public let eta: Vector<3, Scalar>

        @inlinable
        public init(xi: Vector<3, Scalar>, eta: Vector<3, Scalar>) {
            self.xi = xi
            self.eta = eta
        }
    }
}

extension Galerkin.Jacobian.Surface {
    /// The Jacobian of the face with these corners, in the reference cell's order, where `shape` was evaluated.
    @inlinable
    public init(_ shape: Galerkin.Shape<Scalar>, corners: [Point<3, Scalar>]) {
        precondition(corners.count == shape.count, "one corner per shape function")
        var xi = Vector<3, Scalar>.zero
        var eta = Vector<3, Scalar>.zero
        for index in corners.indices {
            let corner = corners[index].coordinates
            xi += corner.scaled(by: shape.xi[index])
            eta += corner.scaled(by: shape.eta[index])
        }
        self.init(xi: xi, eta: eta)
    }

    /// x,ξ × x,η: normal to the face, pointing to the side from which the reference cell's corners turn
    /// counterclockwise, as long as the ratio of the face's area to the reference cell's at this point.
    @inlinable
    public var normal: Vector<3, Scalar> { xi.cross(eta) }

    /// |x,ξ × x,η|: the ratio of the face's area to the reference cell's at this point, never negative, as a surface
    /// has no orientation of its own.
    @inlinable
    public var determinant: Scalar { normal.length }
}

extension Galerkin.Jacobian.Surface: Equatable {}
extension Galerkin.Jacobian.Surface: Hashable where Scalar: Hashable {}
extension Galerkin.Jacobian.Surface: Sendable where Scalar: Sendable {}
