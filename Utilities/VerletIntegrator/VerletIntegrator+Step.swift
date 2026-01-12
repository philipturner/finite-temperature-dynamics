extension VerletIntegrator {
  typealias Singlepoint = (_ positions: [SIMD3<Float>]) -> (
    forces: [SIMD3<Float>],
    energy: Double
  )
  
  // Perform one integration step.
  mutating func step(
    time: Float,
    singlepoint: Singlepoint
  ) {
    if forces == nil {
      forces = singlepoint(positions).forces
    }
    
    integrateVelocities(time: time / 2)
    integratePositions(time: time)
    
    let singlepointResult = singlepoint(positions)
    forces = singlepointResult.forces
    integrateVelocities(time: time / 2)
    
    var maximumForce: Float = .zero
    for atomID in masses.indices {
      let mass = masses[atomID]
      guard mass > 0 else {
        continue
      }
      
      let force = forces![atomID]
      let forceMagnitude = (force * force).sum().squareRoot()
      maximumForce = max(maximumForce, forceMagnitude)
    }
    
    let temperature = VerletIntegrator.temperature(
      velocities: velocities,
      masses: masses)
    
    if let timeForConsole {
      print("time: \(Format.timePs(timeForConsole))", terminator: " | ")
    }
    print("energy: \(Format.energy(singlepointResult.energy))", terminator: " | ")
    print("temp: \(Format.temperature(temperature))", terminator: " | ")
    print("max force: \(Format.force(maximumForce))", terminator: " | ")
    print("Δt: \(Format.time(time))", terminator: " | ")
    print()
  }
  
  mutating func integrateVelocities(time: Float) {
    guard let forces else {
      fatalError("Forces were not generated yet.")
    }
    
    for atomID in masses.indices {
      let mass = masses[atomID]
      guard mass > 0 else {
        continue
      }
      
      let force = forces[atomID]
      velocities[atomID] += time * force / mass
    }
  }
  
  mutating func integratePositions(time: Float) {
    for atomID in masses.indices {
      let velocity = velocities[atomID]
      positions[atomID] += time * velocity
    }
  }
}
