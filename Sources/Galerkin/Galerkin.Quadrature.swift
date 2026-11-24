internal import Numeric

extension Galerkin {
    /// Points of an `N`-dimensional reference cell and their weights: Σ wᵢ f(pᵢ) is the integral of f over the
    /// cell, exactly when f is a polynomial of total degree at most `degree`.
    public struct Quadrature<let N: Int, Scalar: BinaryFloatingPoint> {
        /// The points, in the reference cell's coordinates (ξ, η) or (ξ, η, ζ).
        public let points: [Point<N, Scalar>]
        /// The weight of each point.
        public let weights: [Scalar]
        /// The highest total degree integrated exactly.
        public let degree: Int

        @inlinable
        public init(points: [Point<N, Scalar>], weights: [Scalar], degree: Int) {
            precondition(points.count == weights.count, "one weight per point")
            self.points = points
            self.weights = weights
            self.degree = degree
        }
    }
}

extension Galerkin.Quadrature {
    /// Σ wᵢ f(pᵢ).
    @inlinable
    public func integral(_ integrand: (Point<N, Scalar>) -> Scalar) -> Scalar {
        var sum: Scalar = 0
        for (point, weight) in zip(points, weights) { sum += weight * integrand(point) }
        return sum
    }
}

extension Galerkin.Quadrature where N == 1 {
    /// The `count`-point Gauss–Legendre rule on [−1, 1], exact for polynomials of degree 2 `count` − 1, its nodes
    /// ascending.
    ///
    /// Each node is the root of the Legendre polynomial Pₙ found by Newton's method from Tricomi's estimate
    /// cos(π (i + ¾) / (n + ½)), and its weight is 2 / ((1 − x²) Pₙ′(x)²).
    public static func gauss(count: Int) -> (nodes: [Scalar], weights: [Scalar]) {
        precondition(count > 0, "a Gauss rule has at least one point")
        var nodes = [Scalar](repeating: 0, count: count)
        var weights = [Scalar](repeating: 0, count: count)
        /// Pₙ(x) and Pₙ′(x), from P₀ = 1, P₁ = x and (k + 1) Pₖ₊₁ = (2k + 1) x Pₖ − k Pₖ₋₁.
        func legendre(_ x: Scalar) -> (value: Scalar, derivative: Scalar) {
            var previous: Scalar = 1
            var current = x
            var k = 1
            while k < count {
                let next = (Scalar(2 * k + 1) * x * current - Scalar(k) * previous) / Scalar(k + 1)
                previous = current
                current = next
                k += 1
            }
            return (current, Scalar(count) * (x * current - previous) / (x * x - 1))
        }
        for index in 0..<(count + 1) / 2 {
            var x = Scalar(Double.cos(Double.pi * (Double(index) + 0.75) / (Double(count) + 0.5)))
            var step = 0
            while step < 100 {
                let p = legendre(x)
                let moved = x - p.value / p.derivative
                step += 1
                let settled = Swift.abs(moved - x) <= 2 * Scalar.ulpOfOne
                x = moved
                if settled { break }
            }
            let derivative = legendre(x).derivative
            let weight = 2 / ((1 - x * x) * derivative * derivative)
            nodes[index] = -x
            nodes[count - 1 - index] = x
            weights[index] = weight
            weights[count - 1 - index] = weight
        }
        if count % 2 == 1 { nodes[count / 2] = 0 }
        return (nodes, weights)
    }
}

extension Galerkin.Quadrature: Equatable {}
extension Galerkin.Quadrature: Hashable where Scalar: Hashable {}
extension Galerkin.Quadrature: Sendable where Scalar: Sendable {}
