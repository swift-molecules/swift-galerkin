import Galerkin
import Testing

extension `Galerkin Tests` {
    @Suite struct Surface {}
}

extension `Galerkin Tests`.Surface {
    static let flat = [Point<2, Double>(x: 0, y: 0), Point(x: 2, y: 0.3), Point(x: 1.7, y: 1.6), Point(x: -0.2, y: 1.1)]

    /// `flat` turned out of every coordinate plane, about x by the angle whose cosine is 3/5 and about y by the one
    /// whose cosine is 5/13, and lifted.
    static let tilted: [Point<3, Double>] = {
        let (ca, sa, cb, sb) = (0.6, 0.8, 5.0 / 13, 12.0 / 13)
        return flat.map { p in
            // About x: (x, y cos a, y sin a); then about y: (x cos b + z sin b, y, −x sin b + z cos b).
            let (x, y, z) = (p.x, p.y * ca, p.y * sa)
            return Point(x: x * cb + z * sb + 0.3, y: y - 0.2, z: -x * sb + z * cb + 1)
        }
    }()

    @Test
    func `A flat face's area does not change when it is turned in space`() {
        let rule = Galerkin.Reference.Quadrilateral.quadrature(degree: 2, Double.self)
        let area = rule.integral {
            Galerkin.Jacobian.Surface(Galerkin.Reference.Quadrilateral.shape(at: $0), corners: Self.tilted).determinant
        }
        // The shoelace area of the face before it was turned.
        var shoelace = 0.0
        for i in 0..<4 {
            let (a, b) = (Self.flat[i], Self.flat[(i + 1) % 4])
            shoelace += a.x * b.y - b.x * a.y
        }
        #expect(abs(area - shoelace / 2) < 1e-13, "\(area)")
        let triangle = [Point<3, Double>(x: 1, y: 1, z: 0), Point(x: 4, y: 2, z: 0), Point(x: 2, y: 5, z: 0)]
        let found = Galerkin.Reference.Triangle.quadrature(degree: 0, Double.self).integral {
            Galerkin.Jacobian.Surface(Galerkin.Reference.Triangle.shape(at: $0), corners: triangle).determinant
        }
        #expect(abs(found - 5.5) < 1e-14)
    }

    @Test
    func `The normal is square to the face and points where its corners turn counterclockwise`() {
        let rule = Galerkin.Reference.Quadrilateral.quadrature(degree: 3, Double.self)
        let shapes = rule.points.map { Galerkin.Reference.Quadrilateral.shape(at: $0) }
        let first = Galerkin.Jacobian.Surface(shapes[0], corners: Self.tilted).normal
        for shape in shapes {
            let jacobian = Galerkin.Jacobian.Surface(shape, corners: Self.tilted)
            let normal = jacobian.normal
            #expect(abs(normal.dot(jacobian.xi)) < 1e-14)
            #expect(abs(normal.dot(jacobian.eta)) < 1e-14)
            // A planar face has one normal direction everywhere.
            #expect(abs(normal.dot(first) / (normal.length * first.length) - 1) < 1e-14)
        }
        let square = [
            Point<3, Double>(x: 0, y: 0, z: 2), Point(x: 1, y: 0, z: 2), Point(x: 1, y: 1, z: 2), Point(x: 0, y: 1, z: 2),
        ]
        let up = Galerkin.Jacobian.Surface(Galerkin.Reference.Quadrilateral.shape(at: Point(x: 0.3, y: -0.2)), corners: square)
        #expect(up.normal.z > 0 && abs(up.normal.x) < 1e-15 && abs(up.normal.y) < 1e-15)
        let down = Galerkin.Jacobian.Surface(
            Galerkin.Reference.Quadrilateral.shape(at: Point(x: 0.3, y: -0.2)),
            corners: square.reversed()
        )
        #expect(down.normal.z < 0)
        #expect(abs(down.determinant - up.determinant) < 1e-15)
    }

    @Test
    func `A warped face's area is the area of the surface its corners span`() {
        // The corners (0, 0, 0), (1, 0, 0), (1, 1, 1), (0, 1, 0) span the saddle z = x y over the unit square, which
        // the bilinear map gives exactly. Its area is ∫∫ √(1 + x² + y²) dx dy, here by Simpson's rule on a fine grid.
        let corners = [
            Point<3, Double>(x: 0, y: 0, z: 0), Point(x: 1, y: 0, z: 0), Point(x: 1, y: 1, z: 1), Point(x: 0, y: 1, z: 0),
        ]
        let steps = 400
        func weight(_ i: Int) -> Double { i == 0 || i == steps ? 1 : (i % 2 == 1 ? 4 : 2) }
        var simpson = 0.0
        for i in 0...steps {
            for j in 0...steps {
                let (x, y) = (Double(i) / Double(steps), Double(j) / Double(steps))
                simpson += weight(i) * weight(j) * (1 + x * x + y * y).squareRoot()
            }
        }
        let exact = simpson / Double(9 * steps * steps)
        let rule = Galerkin.Reference.Quadrilateral.quadrature(degree: 13, Double.self)
        let area = rule.integral {
            Galerkin.Jacobian.Surface(Galerkin.Reference.Quadrilateral.shape(at: $0), corners: corners).determinant
        }
        #expect(abs(area - exact) < 1e-9, "\(area) against \(exact)")
    }
}
