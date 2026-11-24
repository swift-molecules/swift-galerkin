// swift-tools-version: 6.4

import PackageDescription

let package = Package(
    name: "swift-galerkin",
    platforms: [
        .macOS(.v27),
        .iOS(.v27),
        .tvOS(.v27),
        .watchOS(.v27),
        .visionOS(.v27),
    ],
    products: [
        .library(name: "Galerkin", targets: ["Galerkin"])
    ],
    dependencies: [
        .package(url: "https://github.com/swift-molecules/swift-sparse.git", branch: "main"),
        .package(url: "https://github.com/swift-atoms/swift-point.git", branch: "main"),
        .package(url: "https://github.com/swift-atoms/swift-numeric.git", branch: "main"),
    ],
    targets: [
        .target(
            name: "Galerkin",
            dependencies: [
                .product(name: "Sparse", package: "swift-sparse"),
                .product(name: "Point", package: "swift-point"),
                .product(name: "Numeric", package: "swift-numeric"),
            ]
        ),
        .testTarget(
            name: "Galerkin Tests",
            dependencies: [
                .target(name: "Galerkin")
            ]
        ),
    ],
    swiftLanguageModes: [.v6]
)

for target in package.targets where ![.system, .binary, .plugin, .macro].contains(target.type) {
    target.swiftSettings =
        (target.swiftSettings ?? []) + [
            .strictMemorySafety(),
            .enableUpcomingFeature("ExistentialAny"),
            .enableUpcomingFeature("InternalImportsByDefault"),
            .enableUpcomingFeature("MemberImportVisibility"),
            .enableUpcomingFeature("NonisolatedNonsendingByDefault"),
            .enableUpcomingFeature("InferIsolatedConformances"),
            .enableExperimentalFeature("Lifetimes"),
            .treatAllWarnings(as: .error),
        ]
}
