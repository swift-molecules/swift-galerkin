extension Galerkin.Reference {
    /// The cube [−1, 1]³ with trilinear shape functions: its corners are the bottom face ζ = −1 counterclockwise
    /// from (−1, −1, −1), seen from above, then the top face ζ = 1 in the same order.
    public enum Hexahedron {}
}

extension Galerkin.Reference.Hexahedron {
    /// The corners: (−1, −1, −1), (1, −1, −1), (1, 1, −1), (−1, 1, −1), then the same with ζ = 1.
    @inlinable
    public static func corners<Scalar: BinaryFloatingPoint>(_: Scalar.Type = Scalar.self) -> [Point<3, Scalar>] {
        [
            Point(x: -1, y: -1, z: -1), Point(x: 1, y: -1, z: -1), Point(x: 1, y: 1, z: -1), Point(x: -1, y: 1, z: -1),
            Point(x: -1, y: -1, z: 1), Point(x: 1, y: -1, z: 1), Point(x: 1, y: 1, z: 1), Point(x: -1, y: 1, z: 1),
        ]
    }

    /// The trilinear shape functions Nᵢ = ⅛ (1 + ξᵢ ξ)(1 + ηᵢ η)(1 + ζᵢ ζ) and their derivatives at
    /// `point` = (ξ, η, ζ).
    @inlinable
    public static func shape<Scalar: BinaryFloatingPoint>(at point: Point<3, Scalar>) -> Galerkin.Shape<Scalar> {
        let xi = point.x
        let eta = point.y
        let zeta = point.z
        var values: [Scalar] = []
        var alongXi: [Scalar] = []
        var alongEta: [Scalar] = []
        var alongZeta: [Scalar] = []
        values.reserveCapacity(8)
        alongXi.reserveCapacity(8)
        alongEta.reserveCapacity(8)
        alongZeta.reserveCapacity(8)
        for corner in corners(Scalar.self) {
            let a = 1 + corner.x * xi
            let b = 1 + corner.y * eta
            let c = 1 + corner.z * zeta
            values.append(a * b * c / 8)
            alongXi.append(corner.x * b * c / 8)
            alongEta.append(corner.y * a * c / 8)
            alongZeta.append(corner.z * a * b / 8)
        }
        return Galerkin.Shape(values: values, xi: alongXi, eta: alongEta, zeta: alongZeta)
    }

    /// The tensor-product Gauss rule exact for every polynomial of degree `degree` in each of ξ, η and ζ, and so
    /// for every polynomial of total degree `degree`.
    @inlinable
    public static func quadrature<Scalar: BinaryFloatingPoint>(
        degree: Int,
        _: Scalar.Type = Scalar.self
    ) -> Galerkin.Quadrature<3, Scalar> {
        precondition(degree >= 0, "a degree is not negative")
        let line = Galerkin.Quadrature<1, Scalar>.gauss(count: degree / 2 + 1)
        var points: [Point<3, Scalar>] = []
        var weights: [Scalar] = []
        for (zeta, outer) in zip(line.nodes, line.weights) {
            for (eta, middle) in zip(line.nodes, line.weights) {
                for (xi, inner) in zip(line.nodes, line.weights) {
                    points.append(Point(x: xi, y: eta, z: zeta))
                    weights.append(outer * middle * inner)
                }
            }
        }
        return Galerkin.Quadrature(points: points, weights: weights, degree: degree)
    }
}
