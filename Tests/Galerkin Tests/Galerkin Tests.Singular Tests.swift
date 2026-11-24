import Galerkin
import Numeric
import Testing

extension `Galerkin Tests` {
    @Suite struct Singular {}
}

extension `Galerkin Tests`.Singular {
    /// ∫ 1 / |x − p| dA over the triangle (0, 0), (1, 0), (0, 1) from a point p in it: each edge at distance h from
    /// p, with its ends at signed distances s₁ < s₂ along it from the foot of the perpendicular, gives
    /// h (asinh(s₂ / h) − asinh(s₁ / h)).
    static func inverseDistance(from p: Point<2, Double>) -> Double {
        let corners = Galerkin.Reference.Triangle.corners(Double.self)
        var sum = 0.0
        for edge in 0..<3 {
            let a = corners[edge]
            let b = corners[(edge + 1) % 3]
            let length = ((b.x - a.x) * (b.x - a.x) + (b.y - a.y) * (b.y - a.y)).squareRoot()
            let tx = (b.x - a.x) / length
            let ty = (b.y - a.y) / length
            let h = (b.x - a.x) * (a.y - p.y) / length - (b.y - a.y) * (a.x - p.x) / length
            guard abs(h) > 1e-15 else { continue }
            let s1 = (a.x - p.x) * tx + (a.y - p.y) * ty
            let s2 = (b.x - p.x) * tx + (b.y - p.y) * ty
            sum += abs(h) * (asinh(s2 / abs(h)) - asinh(s1 / abs(h)))
        }
        return sum
    }

    static func asinh(_ x: Double) -> Double {
        x >= 0 ? Double.log(x + (x * x + 1).squareRoot()) : -Double.log(-x + (x * x + 1).squareRoot())
    }

    static let singular: [Point<2, Double>] = [
        Point(x: 1.0 / 3, y: 1.0 / 3), Point(x: 0.1, y: 0.7), Point(x: 0.5, y: 0), Point(x: 0, y: 0),
        Point(x: 0.2, y: 0.8),
    ]

    @Test(arguments: 0...8)
    func `The rule around a point integrates every polynomial of its degree exactly`(degree: Int) {
        for point in Self.singular {
            let rule = Galerkin.Reference.Triangle.quadrature(degree: degree, around: point)
            #expect(rule.degree == degree)
            for a in 0...degree {
                for b in 0...(degree - a) {
                    let found = rule.integral { `Galerkin Tests`.power($0.x, a) * `Galerkin Tests`.power($0.y, b) }
                    let exact = `Galerkin Tests`.triangle(a, b)
                    #expect(abs(found - exact) < 1e-14, "around \(point), ξ^\(a) η^\(b): \(found), exact \(exact)")
                }
            }
        }
    }

    @Test
    func `The rule around a point converges on 1 over the distance to it`() {
        for point in Self.singular {
            let exact = Self.inverseDistance(from: point)
            var errors: [Double] = []
            for degree in [4, 12, 24] {
                let rule = Galerkin.Reference.Triangle.quadrature(degree: degree, around: point)
                let found = rule.integral { x in
                    1 / ((x.x - point.x) * (x.x - point.x) + (x.y - point.y) * (x.y - point.y)).squareRoot()
                }
                errors.append(abs(found - exact) / exact)
            }
            #expect(errors[1] < errors[0] && errors[2] < errors[1], "around \(point): \(errors)")
            #expect(errors[2] < 1e-3, "around \(point): \(errors)")
        }
    }

    @Test
    func `The rule around a point keeps its points in the triangle and its weights positive`() {
        for point in Self.singular {
            let rule = Galerkin.Reference.Triangle.quadrature(degree: 6, around: point)
            for (x, weight) in zip(rule.points, rule.weights) {
                #expect(x.x >= 0 && x.y >= 0 && x.x + x.y <= 1 + 1e-15)
                #expect(weight > 0)
            }
        }
    }
}
