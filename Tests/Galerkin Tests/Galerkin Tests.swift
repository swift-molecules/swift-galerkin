import Galerkin
import Testing

@Suite struct `Galerkin Tests` {}

extension `Galerkin Tests` {
    /// A reproducible stream of numbers in [−1, 1).
    struct Stream {
        var state: UInt64

        mutating func next() -> Double {
            state = state &* 6_364_136_223_846_793_005 &+ 1_442_695_040_888_963_407
            return Double(state >> 11) / Double(1 << 52) - 1
        }
    }

    /// n!
    static func factorial(_ n: Int) -> Double { n < 2 ? 1 : Double(n) * factorial(n - 1) }

    /// ∫ ξᵃ ηᵇ over [−1, 1]².
    static func square(_ a: Int, _ b: Int) -> Double {
        func line(_ p: Int) -> Double { p % 2 == 1 ? 0 : 2 / Double(p + 1) }
        return line(a) * line(b)
    }

    /// ∫ ξᵃ ηᵇ over the triangle (0, 0), (1, 0), (0, 1): a! b! / (a + b + 2)!.
    static func triangle(_ a: Int, _ b: Int) -> Double {
        factorial(a) * factorial(b) / factorial(a + b + 2)
    }

    static func power(_ x: Double, _ p: Int) -> Double {
        var result = 1.0
        for _ in 0..<p { result *= x }
        return result
    }
}

extension `Galerkin Tests` {
    @Suite struct Quadrature {}
}

extension `Galerkin Tests`.Quadrature {
    @Test(arguments: 0...11)
    func `The square's rule integrates every polynomial of its degree exactly`(degree: Int) {
        let rule = Galerkin.Reference.Quadrilateral.quadrature(degree: degree, Double.self)
        #expect(rule.degree == degree)
        for a in 0...degree {
            for b in 0...(degree - a) {
                let found = rule.integral { `Galerkin Tests`.power($0.x, a) * `Galerkin Tests`.power($0.y, b) }
                let exact = `Galerkin Tests`.square(a, b)
                #expect(abs(found - exact) < 1e-13, "ξ^\(a) η^\(b): \(found), exact \(exact)")
            }
        }
    }

    @Test(arguments: 0...11)
    func `The triangle's rule integrates every polynomial of its degree exactly`(degree: Int) {
        let rule = Galerkin.Reference.Triangle.quadrature(degree: degree, Double.self)
        #expect(rule.degree == degree)
        for a in 0...degree {
            for b in 0...(degree - a) {
                let found = rule.integral { `Galerkin Tests`.power($0.x, a) * `Galerkin Tests`.power($0.y, b) }
                let exact = `Galerkin Tests`.triangle(a, b)
                #expect(abs(found - exact) < 1e-14, "ξ^\(a) η^\(b): \(found), exact \(exact)")
            }
        }
    }

    @Test(arguments: 0...11)
    func `The triangle's points lie inside it and its weights are positive`(degree: Int) {
        let rule = Galerkin.Reference.Triangle.quadrature(degree: degree, Double.self)
        for (point, weight) in zip(rule.points, rule.weights) {
            #expect(point.x > 0 && point.y > 0 && point.x + point.y < 1)
            #expect(weight > 0)
        }
    }

    @Test
    func `The square's rule is no better than its degree when that is odd`() {
        // n Gauss points integrate ξ^(2n) wrongly, so the rule of degree 2n − 1 misses degree 2n.
        for count in 1...6 {
            let rule = Galerkin.Reference.Quadrilateral.quadrature(degree: 2 * count - 1, Double.self)
            let found = rule.integral { `Galerkin Tests`.power($0.x, 2 * count) }
            #expect(abs(found - `Galerkin Tests`.square(2 * count, 0)) > 1e-6)
        }
    }

    @Test
    func `Gauss nodes and weights match the closed forms`() {
        let two = Galerkin.Quadrature<Double>.gauss(count: 2)
        #expect(abs(two.nodes[1] - 1 / 3.0.squareRoot()) < 1e-15 && two.nodes[0] == -two.nodes[1])
        #expect(two.weights.allSatisfy { abs($0 - 1) < 1e-15 })
        let three = Galerkin.Quadrature<Double>.gauss(count: 3)
        #expect(three.nodes[1] == 0)
        #expect(abs(three.nodes[2] - 0.6.squareRoot()) < 1e-15)
        #expect(abs(three.weights[1] - 8 / 9.0) < 1e-15 && abs(three.weights[2] - 5 / 9.0) < 1e-15)
        let float = Galerkin.Quadrature<Float>.gauss(count: 3)
        #expect(abs(float.nodes[2] - Float(0.6).squareRoot()) < 1e-6)
    }
}

extension `Galerkin Tests` {
    @Suite struct Shape {}
}

extension `Galerkin Tests`.Shape {
    static func points(_ count: Int, inTriangle: Bool) -> [Point<2, Double>] {
        var stream = `Galerkin Tests`.Stream(state: 7)
        return (0..<count).map { _ in
            var x = stream.next()
            var y = stream.next()
            if inTriangle {
                x = (x + 1) / 2
                y = (y + 1) / 2
                if x + y > 1 { (x, y) = (1 - x, 1 - y) }
            }
            return Point(x: x, y: y)
        }
    }

    @Test
    func `The square's shape functions are a partition of unity`() {
        for point in Self.points(64, inTriangle: false) {
            let shape = Galerkin.Reference.Quadrilateral.shape(at: point)
            #expect(abs(shape.values.reduce(0, +) - 1) < 1e-15)
            #expect(abs(shape.xi.reduce(0, +)) < 1e-15)
            #expect(abs(shape.eta.reduce(0, +)) < 1e-15)
        }
    }

    @Test
    func `The triangle's shape functions are a partition of unity`() {
        for point in Self.points(64, inTriangle: true) {
            let shape = Galerkin.Reference.Triangle.shape(at: point)
            #expect(abs(shape.values.reduce(0, +) - 1) < 1e-15)
            #expect(shape.xi.reduce(0, +) == 0)
            #expect(shape.eta.reduce(0, +) == 0)
        }
    }

    @Test
    func `Each shape function is one at its own corner and zero at the others`() {
        let square = Galerkin.Reference.Quadrilateral.corners(Double.self)
        for (index, corner) in square.enumerated() {
            let values = Galerkin.Reference.Quadrilateral.shape(at: corner).values
            #expect(values == values.indices.map { $0 == index ? 1 : 0 })
        }
        let triangle = Galerkin.Reference.Triangle.corners(Double.self)
        for (index, corner) in triangle.enumerated() {
            let values = Galerkin.Reference.Triangle.shape(at: corner).values
            #expect(values == values.indices.map { $0 == index ? 1 : 0 })
        }
    }

    @Test
    func `Shape functions reproduce a linear field and its gradient on a distorted element`() {
        let corners = [
            Point<2, Double>(x: 0, y: 0), Point(x: 2, y: 0.3), Point(x: 1.7, y: 1.6), Point(x: -0.2, y: 1.1),
        ]
        let field = { (p: Point<2, Double>) in 3 + 2 * p.x - 5 * p.y }
        let values = corners.map(field)
        for point in Self.points(16, inTriangle: false) {
            let shape = Galerkin.Reference.Quadrilateral.shape(at: point)
            let jacobian = Galerkin.Jacobian(shape, corners: corners)
            #expect(jacobian.determinant > 0)
            let mapped = Point<2, Double>(
                x: shape.interpolated(corners.map(\.x)),
                y: shape.interpolated(corners.map(\.y))
            )
            #expect(abs(shape.interpolated(values) - field(mapped)) < 1e-13)
            let gradients = jacobian.gradients(shape)
            let gx = zip(gradients.x, values).reduce(0) { $0 + $1.0 * $1.1 }
            let gy = zip(gradients.y, values).reduce(0) { $0 + $1.0 * $1.1 }
            #expect(abs(gx - 2) < 1e-12 && abs(gy + 5) < 1e-12)
        }
    }

    @Test
    func `The determinant integrated over the reference cell is the element's area`() {
        let quad = [Point<2, Double>(x: 0, y: 0), Point(x: 2, y: 0.3), Point(x: 1.7, y: 1.6), Point(x: -0.2, y: 1.1)]
        // The shoelace area.
        var shoelace = 0.0
        for i in 0..<4 { shoelace += quad[i].x * quad[(i + 1) % 4].y - quad[(i + 1) % 4].x * quad[i].y }
        let rule = Galerkin.Reference.Quadrilateral.quadrature(degree: 2, Double.self)
        let area = rule.integral {
            Galerkin.Jacobian(Galerkin.Reference.Quadrilateral.shape(at: $0), corners: quad).determinant
        }
        #expect(abs(area - shoelace / 2) < 1e-14)
        let triangle = [Point<2, Double>(x: 1, y: 1), Point(x: 4, y: 2), Point(x: 2, y: 5)]
        let found = Galerkin.Reference.Triangle.quadrature(degree: 0, Double.self).integral {
            Galerkin.Jacobian(Galerkin.Reference.Triangle.shape(at: $0), corners: triangle).determinant
        }
        #expect(abs(found - 5.5) < 1e-14)
    }
}

extension `Galerkin Tests` {
    @Suite struct Assembly {}
}

extension `Galerkin Tests`.Assembly {
    @Test
    func `Unknowns are numbered vertex by vertex, skipping the held ones`() {
        let unknowns = Galerkin.Unknowns(vertices: 3, components: 2, held: [0: [1], 2: [0, 1]])
        #expect(unknowns.count == 3)
        #expect(unknowns.vertices == 3)
        #expect(unknowns[0, 0] == 0 && unknowns[0, 1] == nil)
        #expect(unknowns[1, 0] == 1 && unknowns[1, 1] == 2)
        #expect(unknowns[2, 0] == nil && unknowns[2, 1] == nil)
    }

    /// A random symmetric `order` × `order` matrix, row-major.
    static func symmetric(_ order: Int, _ stream: inout `Galerkin Tests`.Stream) -> [Double] {
        var matrix = [Double](repeating: 0, count: order * order)
        for a in 0..<order {
            for b in a..<order {
                let value = stream.next()
                matrix[order * a + b] = value
                matrix[order * b + a] = value
            }
        }
        return matrix
    }

    @Test
    func `Assembly equals the dense sum of the turned element matrices over the free unknowns`() {
        // A 3 × 2 grid of quads with three components at each of its 12 vertices.
        let columns = 3
        let components = 3
        var elements: [[Int]] = []
        for row in 0..<2 {
            for column in 0..<columns {
                let first = row * (columns + 1) + column
                elements.append([first, first + 1, first + columns + 2, first + columns + 1])
            }
        }
        let vertices = 3 * (columns + 1)
        let held: [Int: Set<Int>] = [0: [0, 1, 2], 1: [0], 5: [0, 2], 11: [1]]
        let rotations: [Int: Galerkin.Rotation<Double>] = [
            1: Galerkin.Rotation(1, 2, along: (x: 3, y: 4)),
            6: Galerkin.Rotation(0, 2, along: (x: -1, y: 2)),
            11: Galerkin.Rotation(1, 2, along: (x: 0.6, y: -0.8)),
        ]
        let unknowns = Galerkin.Unknowns(vertices: vertices, components: components, held: held)
        var assembly = Galerkin.Assembly(unknowns, rotations: rotations)
        var stream = `Galerkin Tests`.Stream(state: 42)
        let matrices = elements.map { _ in Self.symmetric(4 * components, &stream) }
        for (corners, matrix) in zip(elements, matrices) { assembly.add(matrix, at: corners) }
        let sparse = assembly.matrix

        // The dense sum over every component, then Tᵀ A T with T the block-diagonal turn, then the free rows.
        let full = vertices * components
        var dense = [Double](repeating: 0, count: full * full)
        for (corners, matrix) in zip(elements, matrices) {
            for (i, p) in corners.enumerated() {
                for (j, q) in corners.enumerated() {
                    for a in 0..<components {
                        for b in 0..<components {
                            dense[full * (components * p + a) + components * q + b] +=
                                matrix[4 * components * (components * i + a) + components * j + b]
                        }
                    }
                }
            }
        }
        var turn = [Double](repeating: 0, count: full * full)
        for index in 0..<full { turn[full * index + index] = 1 }
        for (vertex, rotation) in rotations {
            let p = components * vertex + rotation.first
            let q = components * vertex + rotation.second
            turn[full * p + p] = rotation.cosine
            turn[full * p + q] = -rotation.sine
            turn[full * q + p] = rotation.sine
            turn[full * q + q] = rotation.cosine
        }
        var turned = [Double](repeating: 0, count: full * full)
        for r in 0..<full {
            for c in 0..<full {
                var sum = 0.0
                for k in 0..<full {
                    for l in 0..<full where turn[full * k + r] != 0 && turn[full * l + c] != 0 {
                        sum += turn[full * k + r] * dense[full * k + l] * turn[full * l + c]
                    }
                }
                turned[full * r + c] = sum
            }
        }
        var free: [Int] = []
        for vertex in 0..<vertices {
            for component in 0..<components where !(held[vertex]?.contains(component) ?? false) {
                free.append(components * vertex + component)
            }
        }
        #expect(sparse.order == free.count)
        #expect(sparse.order == unknowns.count)
        for (r, row) in free.enumerated() {
            for (c, column) in free.enumerated() {
                let found = sparse[Ordinal(UInt(r)), Ordinal(UInt(c))]
                #expect(abs(found - turned[full * row + column]) < 1e-13, "\(r), \(c)")
            }
        }
    }

    @Test
    func `Turning a component by the identity frame changes nothing`() {
        var stream = `Galerkin Tests`.Stream(state: 3)
        let matrix = Self.symmetric(6, &stream)
        var turned = matrix
        Galerkin.Rotation(0, 1, along: (x: 2, y: 0)).turn(&turned, order: 6, offset: 3)
        #expect(turned == matrix)
    }
}
