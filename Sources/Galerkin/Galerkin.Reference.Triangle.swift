extension Galerkin.Reference {
    /// The triangle with corners (0, 0), (1, 0), (0, 1), in that order, with linear shape functions.
    public enum Triangle {}
}

extension Galerkin.Reference.Triangle {
    /// The corners (0, 0), (1, 0), (0, 1).
    @inlinable
    public static func corners<Scalar: BinaryFloatingPoint>(_: Scalar.Type = Scalar.self) -> [Point<2, Scalar>] {
        [Point(x: 0, y: 0), Point(x: 1, y: 0), Point(x: 0, y: 1)]
    }

    /// The linear shape functions 1 − ξ − η, ξ and η, and their derivatives, at `point` = (ξ, η).
    @inlinable
    public static func shape<Scalar: BinaryFloatingPoint>(at point: Point<2, Scalar>) -> Galerkin.Shape<Scalar> {
        Galerkin.Shape(
            values: [1 - point.x - point.y, point.x, point.y],
            xi: [-1, 1, 0],
            eta: [-1, 0, 1]
        )
    }

    /// A Gauss rule exact for every polynomial of total degree `degree`.
    ///
    /// It is the square's Gauss rule collapsed onto the triangle (Duffy): ξ = s (1 − t), η = t, dξ dη = (1 − t) ds dt,
    /// with one more point along t to carry the factor 1 − t. Every point lies inside the triangle and every weight is
    /// positive.
    @inlinable
    public static func quadrature<Scalar: BinaryFloatingPoint>(
        degree: Int,
        _: Scalar.Type = Scalar.self
    ) -> Galerkin.Quadrature<2, Scalar> {
        precondition(degree >= 0, "a degree is not negative")
        let along = Galerkin.Quadrature<1, Scalar>.gauss(count: degree / 2 + 1)
        let across = Galerkin.Quadrature<1, Scalar>.gauss(count: (degree + 1) / 2 + 1)
        var points: [Point<2, Scalar>] = []
        var weights: [Scalar] = []
        for (v, weight) in zip(across.nodes, across.weights) {
            let t = (1 + v) / 2
            for (u, other) in zip(along.nodes, along.weights) {
                let s = (1 + u) / 2
                points.append(Point(x: s * (1 - t), y: t))
                weights.append(weight * other * (1 - t) / 4)
            }
        }
        return Galerkin.Quadrature(points: points, weights: weights, degree: degree)
    }
}
