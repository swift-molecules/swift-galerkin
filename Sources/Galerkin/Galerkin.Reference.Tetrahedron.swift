extension Galerkin.Reference {
    /// The tetrahedron with corners (0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), in that order, with linear shape
    /// functions.
    public enum Tetrahedron {}
}

extension Galerkin.Reference.Tetrahedron {
    /// The corners (0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1).
    @inlinable
    public static func corners<Scalar: BinaryFloatingPoint>(_: Scalar.Type = Scalar.self) -> [Point<3, Scalar>] {
        [Point(x: 0, y: 0, z: 0), Point(x: 1, y: 0, z: 0), Point(x: 0, y: 1, z: 0), Point(x: 0, y: 0, z: 1)]
    }

    /// The linear shape functions 1 − ξ − η − ζ, ξ, η and ζ, and their derivatives, at `point` = (ξ, η, ζ).
    @inlinable
    public static func shape<Scalar: BinaryFloatingPoint>(at point: Point<3, Scalar>) -> Galerkin.Shape<Scalar> {
        Galerkin.Shape(
            values: [1 - point.x - point.y - point.z, point.x, point.y, point.z],
            xi: [-1, 1, 0, 0],
            eta: [-1, 0, 1, 0],
            zeta: [-1, 0, 0, 1]
        )
    }

    /// A Gauss rule exact for every polynomial of total degree `degree`.
    ///
    /// It is the cube's Gauss rule collapsed onto the tetrahedron (Duffy): ξ = a (1 − b)(1 − c), η = b (1 − c),
    /// ζ = c for (a, b, c) in [0, 1]³, with dξ dη dζ = (1 − b)(1 − c)² da db dc, and one more point along b and two
    /// more along c to carry that factor. Every point lies inside the tetrahedron and every weight is positive.
    @inlinable
    public static func quadrature<Scalar: BinaryFloatingPoint>(
        degree: Int,
        _: Scalar.Type = Scalar.self
    ) -> Galerkin.Quadrature<3, Scalar> {
        precondition(degree >= 0, "a degree is not negative")
        let first = Galerkin.Quadrature<1, Scalar>.gauss(count: degree / 2 + 1)
        let second = Galerkin.Quadrature<1, Scalar>.gauss(count: (degree + 1) / 2 + 1)
        let third = Galerkin.Quadrature<1, Scalar>.gauss(count: (degree + 2) / 2 + 1)
        var points: [Point<3, Scalar>] = []
        var weights: [Scalar] = []
        for (w, outer) in zip(third.nodes, third.weights) {
            let c = (1 + w) / 2
            for (v, middle) in zip(second.nodes, second.weights) {
                let b = (1 + v) / 2
                for (u, inner) in zip(first.nodes, first.weights) {
                    let a = (1 + u) / 2
                    points.append(Point(x: a * (1 - b) * (1 - c), y: b * (1 - c), z: c))
                    weights.append(outer * middle * inner * (1 - b) * (1 - c) * (1 - c) / 8)
                }
            }
        }
        return Galerkin.Quadrature(points: points, weights: weights, degree: degree)
    }
}
