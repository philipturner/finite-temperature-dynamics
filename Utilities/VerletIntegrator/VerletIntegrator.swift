struct VerletIntegratorDescriptor {
  /// Required. The mass of each atom (in yoctograms).
  ///
  /// Anchors are specified by setting mass to zero.
  var masses: [Float]?
  
  /// Required. The largest time step that may be taken during simulation, in
  /// picoseconds. Some steps may have a smaller duration.
  var timeStep: Double = 0.0025
}

struct VerletIntegrator {
  // Integrator settings.
  let masses: [Float]
  let timeStep: Double
  
  // Dynamical variables.
  var forces: [SIMD3<Float>]?
  var positions: [SIMD3<Float>]
  var velocities: [SIMD3<Float>]
  
  var timeForConsole: Double?
  
  init(descriptor: VerletIntegratorDescriptor) {
    guard let masses = descriptor.masses else {
      fatalError("Descriptor was incomplete.")
    }
    self.masses = masses
    self.timeStep = descriptor.timeStep
    
    positions = [SIMD3<Float>](repeating: .zero, count: masses.count)
    velocities = [SIMD3<Float>](repeating: .zero, count: masses.count)
  }
  
  /// Simulate the system's evolution for the specified time interval.
  ///
  /// - Parameter time: The time interval, in picoseconds.
  mutating func simulate(
    time: Double,
    singlepoint: Singlepoint
  ) {
    if time == 0 {
      return
    }
    
    // Check whether the arguments are valid.
    guard time > 0, timeStep > 0 else {
      fatalError("Time or time step was invalid.")
    }
    
    // Create rough estimate of step count.
    var stepCountQuotient = (time / timeStep).rounded(.down)
    var stepCountRemainder = (time / timeStep) - stepCountQuotient
    
    guard stepCountQuotient >= 0, stepCountRemainder >= 0 else {
      fatalError("This should never happen.")
    }
    
    // Correct for overshoot and undershoot from floating-point error.
    let epsilon: Double = 1e-4
    if stepCountRemainder < epsilon {
      if stepCountQuotient > 0 {
        stepCountQuotient -= 1
        stepCountRemainder += 1
      }
    } else if stepCountRemainder > 1 + epsilon {
      fatalError("This should never happen.")
    }
    
    if stepCountQuotient == 0 {
      step(
        time: Float(time),
        singlepoint: singlepoint)
    } else {
      let conservativeStepCount = Int(exactly: (time / timeStep).rounded(.up))!
      let conservativeStepSize = time / Double(conservativeStepCount)
      
      for _ in 0..<conservativeStepCount {
        step(
          time: Float(conservativeStepSize),
          singlepoint: singlepoint)
      }
    }
  }
}
