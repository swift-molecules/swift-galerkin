extension Galerkin {
    /// The numbering of a model's unknowns: `components` at every vertex, numbered vertex by vertex and component by
    /// component from zero, skipping those that are held.
    public struct Unknowns {
        /// How many unknowns each vertex has before any is held.
        public let components: Int
        /// How many unknowns are numbered: the order of the assembled matrices.
        public let count: Int
        /// The number of each vertex's each component, or −1 where it is held.
        @usableFromInline let numbers: [Int]

        /// Numbers `components` unknowns at each of `vertices` vertices, except the components `held` names for a
        /// vertex.
        @inlinable
        public init(vertices: Int, components: Int, held: [Int: Set<Int>] = [:]) {
            precondition(vertices >= 0 && components > 0, "a vertex has at least one unknown")
            var numbers = [Int](repeating: -1, count: vertices * components)
            var count = 0
            for vertex in 0..<vertices {
                let fixed = held[vertex] ?? []
                for component in 0..<components where !fixed.contains(component) {
                    numbers[components * vertex + component] = count
                    count += 1
                }
            }
            self.components = components
            self.count = count
            self.numbers = numbers
        }
    }
}

extension Galerkin.Unknowns {
    /// How many vertices are numbered.
    @inlinable
    public var vertices: Int { numbers.count / components }

    /// The number of `component` at `vertex`, or nil where it is held.
    @inlinable
    public subscript(vertex: Int, component: Int) -> Int? {
        let number = numbers[components * vertex + component]
        return number >= 0 ? number : nil
    }
}

extension Galerkin.Unknowns: Hashable {}
extension Galerkin.Unknowns: Sendable {}
