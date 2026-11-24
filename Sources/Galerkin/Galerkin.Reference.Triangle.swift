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

    /// A rule for an integrand that is singular as 1 / r at `point`, a point of the triangle (inside it, on an
    /// edge or at a corner), and still exact for every polynomial of total degree `degree`.
    ///
    /// The triangle is cut at `point` into the triangles it makes with each edge, and `quadrature(degree:)` is
    /// mapped onto each with its collapsed corner at `point`. There the rule's factor 1 − t is proportional to the
    /// distance r, so an integrand f / r, with f smooth, times the map's determinant is bounded and the rule
    /// converges. It converges slowly where `point` is close to an edge next to that edge's length, since 1 / r
    /// along the edge is then nearly singular; a singular part known in closed form is best subtracted first. A
    /// cut of no area, from a point on an edge, is left out.
    public static func quadrature<Scalar: BinaryFloatingPoint>(
        degree: Int,
        around point: Point<2, Scalar>
    ) -> Galerkin.Quadrature<2, Scalar> {
        precondition(
            point.x >= 0 && point.y >= 0 && point.x + point.y <= 1,
            "the singular point lies in the triangle"
        )
        let rule = quadrature(degree: degree, Scalar.self)
        let corners = corners(Scalar.self)
        var points: [Point<2, Scalar>] = []
        var weights: [Scalar] = []
        for edge in 0..<3 {
            let a = corners[edge]
            let b = corners[(edge + 1) % 3]
            // The cut (a, b, point) over the reference triangle: twice its area, the map's determinant.
            let ratio = (b.x - a.x) * (point.y - a.y) - (b.y - a.y) * (point.x - a.x)
            guard ratio > 0 else { continue }
            for (reference, weight) in zip(rule.points, rule.weights) {
                let s = reference.x
                let t = reference.y
                let r = 1 - s - t
                points.append(Point(x: r * a.x + s * b.x + t * point.x, y: r * a.y + s * b.y + t * point.y))
                weights.append(weight * ratio)
            }
        }
        return Galerkin.Quadrature(points: points, weights: weights, degree: degree)
    }
}
