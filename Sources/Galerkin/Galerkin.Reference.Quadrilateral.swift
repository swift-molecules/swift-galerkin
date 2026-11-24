extension Galerkin.Reference {
    /// The square [−1, 1] × [−1, 1] with bilinear shape functions, its corners counterclockwise from (−1, −1).
    public enum Quadrilateral {}
}

extension Galerkin.Reference.Quadrilateral {
    /// The corners, counterclockwise from (−1, −1).
    @inlinable
    public static func corners<Scalar: BinaryFloatingPoint>(_: Scalar.Type = Scalar.self) -> [Point<2, Scalar>] {
        [Point(x: -1, y: -1), Point(x: 1, y: -1), Point(x: 1, y: 1), Point(x: -1, y: 1)]
    }

    /// The bilinear shape functions Nᵢ = ¼ (1 + ξᵢ ξ)(1 + ηᵢ η) and their derivatives at `point` = (ξ, η).
    @inlinable
    public static func shape<Scalar: BinaryFloatingPoint>(at point: Point<2, Scalar>) -> Galerkin.Shape<Scalar> {
        let xi = point.x
        let eta = point.y
        var values: [Scalar] = []
        var alongXi: [Scalar] = []
        var alongEta: [Scalar] = []
        values.reserveCapacity(4)
        alongXi.reserveCapacity(4)
        alongEta.reserveCapacity(4)
        for corner in corners(Scalar.self) {
            let a = corner.x
            let b = corner.y
            values.append((1 + a * xi) * (1 + b * eta) / 4)
            alongXi.append(a * (1 + b * eta) / 4)
            alongEta.append(b * (1 + a * xi) / 4)
        }
        return Galerkin.Shape(values: values, xi: alongXi, eta: alongEta)
    }

    /// The tensor-product Gauss rule exact for every polynomial of degree `degree` in each of ξ and η, and so for
    /// every polynomial of total degree `degree`.
    @inlinable
    public static func quadrature<Scalar: BinaryFloatingPoint>(
        degree: Int,
        _: Scalar.Type = Scalar.self
    ) -> Galerkin.Quadrature<2, Scalar> {
        precondition(degree >= 0, "a degree is not negative")
        let line = Galerkin.Quadrature<1, Scalar>.gauss(count: degree / 2 + 1)
        var points: [Point<2, Scalar>] = []
        var weights: [Scalar] = []
        for (eta, weight) in zip(line.nodes, line.weights) {
            for (xi, other) in zip(line.nodes, line.weights) {
                points.append(Point(x: xi, y: eta))
                weights.append(weight * other)
            }
        }
        return Galerkin.Quadrature(points: points, weights: weights, degree: degree)
    }
}
