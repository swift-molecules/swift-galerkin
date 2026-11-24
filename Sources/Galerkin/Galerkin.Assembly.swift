extension Galerkin {
    /// The sum of element matrices over a model's unknowns, gathered as entries of a `Sparse.Symmetric`.
    ///
    /// Each element matrix is over its corners' components, corner by corner. Adding one turns the components of
    /// every corner that has a rotation, drops the held ones, and adds its upper triangle at the unknowns' numbers.
    public struct Assembly<Scalar: BinaryFloatingPoint> {
        /// The numbering the entries are placed by.
        public let unknowns: Unknowns
        /// The vertices whose components are taken in a turned frame.
        public let rotations: [Int: Rotation<Scalar>]
        /// The entries added so far.
        @usableFromInline var storage: [Sparse.Symmetric<Scalar>.Entry]

        @inlinable
        public init(_ unknowns: Unknowns, rotations: [Int: Rotation<Scalar>] = [:]) {
            self.unknowns = unknowns
            self.rotations = rotations
            self.storage = []
        }
    }
}

extension Galerkin.Assembly {
    /// The entries added so far, each once for a pair of unknowns of one element.
    @inlinable
    public var entries: [Sparse.Symmetric<Scalar>.Entry] { storage }

    /// Adds the symmetric row-major `matrix`, over the components of `corners` in order, corner by corner.
    @inlinable
    public mutating func add(_ matrix: [Scalar], at corners: [Int]) {
        let components = unknowns.components
        let order = components * corners.count
        precondition(matrix.count == order * order, "one row and column per component of each corner")
        var matrix = matrix
        for (local, vertex) in corners.enumerated() {
            rotations[vertex]?.turn(&matrix, order: order, offset: components * local)
        }
        var map: [Int] = []
        map.reserveCapacity(order)
        for vertex in corners {
            for component in 0..<components { map.append(unknowns[vertex, component] ?? -1) }
        }
        for a in 0..<order where map[a] >= 0 {
            for b in a..<order where map[b] >= 0 {
                storage.append(
                    Sparse.Symmetric.Entry(Ordinal(UInt(map[a])), Ordinal(UInt(map[b])), matrix[order * a + b])
                )
            }
        }
    }

    /// The sum of every matrix added, of order `unknowns.count`.
    @inlinable
    public var matrix: Sparse.Symmetric<Scalar> {
        do throws(Sparse.Error) {
            return try Sparse.Symmetric(order: Cardinal(UInt(unknowns.count)), entries: storage)
        } catch {
            preconditionFailure("every unknown is numbered below the order")
        }
    }
}

extension Galerkin.Assembly: Equatable {}
extension Galerkin.Assembly: Sendable where Scalar: Sendable {}
