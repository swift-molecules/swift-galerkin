import Galerkin
import Testing

extension `Galerkin Tests` {
    @Suite struct Solid {}
}

extension `Galerkin Tests`.Solid {
    /// ∫ ξᵃ ηᵇ ζᶜ over [−1, 1]³.
    static func cube(_ a: Int, _ b: Int, _ c: Int) -> Double {
        `Galerkin Tests`.square(a, b) * (c % 2 == 1 ? 0 : 2 / Double(c + 1))
    }

    /// ∫ ξᵃ ηᵇ ζᶜ over the tetrahedron (0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1): a! b! c! / (a + b + c + 3)!.
    static func tetrahedron(_ a: Int, _ b: Int, _ c: Int) -> Double {
        let f = `Galerkin Tests`.factorial
        return f(a) * f(b) * f(c) / f(a + b + c + 3)
    }

    static func monomial(_ p: Point<3, Double>, _ a: Int, _ b: Int, _ c: Int) -> Double {
        let power = `Galerkin Tests`.power
        return power(p.x, a) * power(p.y, b) * power(p.z, c)
    }

    static func points(_ count: Int, inTetrahedron: Bool) -> [Point<3, Double>] {
        var stream = `Galerkin Tests`.Stream(state: 11)
        var result: [Point<3, Double>] = []
        while result.count < count {
            let x = stream.next()
            let y = stream.next()
            let z = stream.next()
            if inTetrahedron {
                let p = Point<3, Double>(x: (x + 1) / 2, y: (y + 1) / 2, z: (z + 1) / 2)
                if p.x + p.y + p.z < 1 { result.append(p) }
            } else {
                result.append(Point(x: x, y: y, z: z))
            }
        }
        return result
    }

    /// A hexahedron with no two faces parallel, its corners in the reference cell's order.
    static let brick: [Point<3, Double>] = [
        Point(x: 0, y: 0, z: 0), Point(x: 2, y: 0.2, z: 0.1), Point(x: 2.3, y: 1.8, z: -0.2),
        Point(x: 0.1, y: 1.5, z: 0),
        Point(x: 0.2, y: -0.1, z: 1.4), Point(x: 1.9, y: 0.1, z: 1.7), Point(x: 2.1, y: 2.0, z: 1.5),
        Point(x: -0.3, y: 1.6, z: 1.3),
    ]

    @Test(arguments: 0...9)
    func `The cube's rule integrates every polynomial of its degree exactly`(degree: Int) {
        let rule = Galerkin.Reference.Hexahedron.quadrature(degree: degree, Double.self)
        #expect(rule.degree == degree)
        for a in 0...degree {
            for b in 0...(degree - a) {
                for c in 0...(degree - a - b) {
                    let found = rule.integral { Self.monomial($0, a, b, c) }
                    #expect(abs(found - Self.cube(a, b, c)) < 1e-13, "ξ^\(a) η^\(b) ζ^\(c): \(found)")
                }
            }
        }
    }

    @Test(arguments: 0...9)
    func `The tetrahedron's rule integrates every polynomial of its degree exactly`(degree: Int) {
        let rule = Galerkin.Reference.Tetrahedron.quadrature(degree: degree, Double.self)
        #expect(rule.degree == degree)
        for a in 0...degree {
            for b in 0...(degree - a) {
                for c in 0...(degree - a - b) {
                    let found = rule.integral { Self.monomial($0, a, b, c) }
                    let exact = Self.tetrahedron(a, b, c)
                    #expect(abs(found - exact) < 1e-15, "ξ^\(a) η^\(b) ζ^\(c): \(found), exact \(exact)")
                }
            }
        }
    }

    @Test(arguments: 0...9)
    func `The tetrahedron's points lie inside it and its weights are positive`(degree: Int) {
        let rule = Galerkin.Reference.Tetrahedron.quadrature(degree: degree, Double.self)
        for (point, weight) in zip(rule.points, rule.weights) {
            #expect(point.x > 0 && point.y > 0 && point.z > 0 && point.x + point.y + point.z < 1)
            #expect(weight > 0)
        }
    }

    @Test
    func `The cube's rule is no better than its degree when that is odd`() {
        for count in 1...5 {
            let rule = Galerkin.Reference.Hexahedron.quadrature(degree: 2 * count - 1, Double.self)
            let found = rule.integral { `Galerkin Tests`.power($0.z, 2 * count) }
            #expect(abs(found - Self.cube(0, 0, 2 * count)) > 1e-6)
        }
    }

    @Test
    func `The solid shape functions are a partition of unity`() {
        for point in Self.points(64, inTetrahedron: false) {
            let shape = Galerkin.Reference.Hexahedron.shape(at: point)
            #expect(abs(shape.values.reduce(0, +) - 1) < 1e-15)
            for derivative in [shape.xi, shape.eta, shape.zeta] { #expect(abs(derivative.reduce(0, +)) < 1e-15) }
        }
        for point in Self.points(64, inTetrahedron: true) {
            let shape = Galerkin.Reference.Tetrahedron.shape(at: point)
            #expect(abs(shape.values.reduce(0, +) - 1) < 1e-15)
            for derivative in [shape.xi, shape.eta, shape.zeta] { #expect(derivative.reduce(0, +) == 0) }
        }
    }

    @Test
    func `Each solid shape function is one at its own corner and zero at the others`() {
        for (index, corner) in Galerkin.Reference.Hexahedron.corners(Double.self).enumerated() {
            let values = Galerkin.Reference.Hexahedron.shape(at: corner).values
            #expect(values == values.indices.map { $0 == index ? 1 : 0 })
        }
        for (index, corner) in Galerkin.Reference.Tetrahedron.corners(Double.self).enumerated() {
            let values = Galerkin.Reference.Tetrahedron.shape(at: corner).values
            #expect(values == values.indices.map { $0 == index ? 1 : 0 })
        }
    }

    @Test
    func `Solid shape functions reproduce a linear field and its gradient on a distorted element`() {
        let field = { (p: Point<3, Double>) in 3 + 2 * p.x - 5 * p.y + 0.5 * p.z }
        let values = Self.brick.map(field)
        for point in Self.points(16, inTetrahedron: false) {
            let shape = Galerkin.Reference.Hexahedron.shape(at: point)
            let jacobian = Galerkin.Jacobian.Solid(shape, corners: Self.brick)
            #expect(jacobian.determinant > 0)
            let mapped = Point<3, Double>(
                x: shape.interpolated(Self.brick.map(\.x)),
                y: shape.interpolated(Self.brick.map(\.y)),
                z: shape.interpolated(Self.brick.map(\.z))
            )
            #expect(abs(shape.interpolated(values) - field(mapped)) < 1e-13)
            let gradients = jacobian.gradients(shape)
            let g = [gradients.x, gradients.y, gradients.z].map { zip($0, values).reduce(0) { $0 + $1.0 * $1.1 } }
            #expect(abs(g[0] - 2) < 1e-12 && abs(g[1] + 5) < 1e-12 && abs(g[2] - 0.5) < 1e-12)
        }
    }

    @Test
    func `The determinant integrated over the reference cell is the element's volume`() {
        let rule = Galerkin.Reference.Hexahedron.quadrature(degree: 3, Double.self)
        // A square frustum, side 2 below and 1 above, height 1.5: h (A₁ + A₂ + √(A₁ A₂)) / 3.
        let frustum = Galerkin.Reference.Hexahedron.corners(Double.self).map {
            let half = $0.z < 0 ? 1.0 : 0.5
            return Point<3, Double>(x: half * $0.x, y: half * $0.y, z: 0.75 * ($0.z + 1))
        }
        let frustumVolume = rule.integral {
            Galerkin.Jacobian.Solid(Galerkin.Reference.Hexahedron.shape(at: $0), corners: frustum).determinant
        }
        #expect(abs(frustumVolume - 1.5 * (4 + 1 + 2) / 3) < 1e-14)
        // A parallelepiped's volume is |a · (b × c)|.
        let box = Galerkin.Reference.Hexahedron.corners(Double.self).map {
            Point<3, Double>(x: 2 * $0.x + 0.5 * $0.y, y: 3 * $0.y + 0.25 * $0.z, z: 1.5 * $0.z)
        }
        let boxVolume = rule.integral {
            Galerkin.Jacobian.Solid(Galerkin.Reference.Hexahedron.shape(at: $0), corners: box).determinant
        }
        #expect(abs(boxVolume - 8 * 2 * 3 * 1.5) < 1e-12)
        // A tetrahedron's volume is a sixth of its edges' triple product.
        let tetrahedron = [
            Point<3, Double>(x: 1, y: 1, z: 1), Point(x: 4, y: 2, z: 1), Point(x: 2, y: 5, z: 1.5),
            Point(x: 1.5, y: 2, z: 4),
        ]
        let found = Galerkin.Reference.Tetrahedron.quadrature(degree: 0, Double.self).integral {
            Galerkin.Jacobian.Solid(Galerkin.Reference.Tetrahedron.shape(at: $0), corners: tetrahedron).determinant
        }
        let e = (1...3).map { tetrahedron[$0].coordinates - tetrahedron[0].coordinates }
        #expect(abs(found - e[0].dot(e[1].cross(e[2])) / 6) < 1e-14)
    }
}
