extension Galerkin {
    /// Two of a vertex's components taken in a turned frame, which is how a constraint along a direction other than
    /// the axes is imposed: hold the turned component, keep the other.
    ///
    /// With the frame's first axis t = (cos, sin) and second n = (−sin, cos), the components (a, b) along the axes are
    /// p t + q n in terms of the turned ones (p, q): a = cos p − sin q, b = sin p + cos q. An element matrix A over the
    /// axes' components becomes Tᵀ A T over the turned ones.
    public struct Rotation<Scalar: BinaryFloatingPoint> {
        /// The components turned: the first becomes the one along t, the second the one along n.
        public let first: Int
        public let second: Int
        /// The frame's first axis, a unit vector.
        public let cosine: Scalar
        public let sine: Scalar

        /// Turns `first` and `second` into the components along `direction` and along its normal; `direction`
        /// need not be a unit vector.
        @inlinable
        public init(_ first: Int, _ second: Int, along direction: (x: Scalar, y: Scalar)) {
            precondition(first != second, "two different components")
            let length = (direction.x * direction.x + direction.y * direction.y).squareRoot()
            precondition(length > 0, "a direction is not zero")
            self.first = first
            self.second = second
            self.cosine = direction.x / length
            self.sine = direction.y / length
        }
    }
}

extension Galerkin.Rotation {
    /// Replaces the `order` × `order` row-major `matrix` by Tᵀ A T, with these components at `offset`.
    @inlinable
    public func turn(_ matrix: inout [Scalar], order: Int, offset: Int) {
        let p = offset + first
        let q = offset + second
        for row in 0..<order {
            let a = matrix[order * row + p]
            let b = matrix[order * row + q]
            matrix[order * row + p] = a * cosine + b * sine
            matrix[order * row + q] = -a * sine + b * cosine
        }
        for column in 0..<order {
            let a = matrix[order * p + column]
            let b = matrix[order * q + column]
            matrix[order * p + column] = a * cosine + b * sine
            matrix[order * q + column] = -a * sine + b * cosine
        }
    }
}

extension Galerkin.Rotation: Equatable {}
extension Galerkin.Rotation: Hashable where Scalar: Hashable {}
extension Galerkin.Rotation: Sendable where Scalar: Sendable {}
