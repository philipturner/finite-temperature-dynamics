import MM4
import OpenMM

extension VerletIntegrator {
  /// Utility function for finding the mass of any element on the periodic
  /// table.
  /// - Returns: Mass, in yoctograms.
  ///
  /// Like the utility from MM4, except hydrogen masses are doubled, and
  /// sulfur masses are zero.
  static func mass(atomicNumber: UInt8) -> Float {
    var output = MM4Parameters.mass(atomicNumber: atomicNumber)
    if atomicNumber == 1 {
      output *= 2
    } else if atomicNumber == 16 {
      output = 0
    }
    return output
  }
  
  /// Create random velocities, according to the Boltzmann distribution.
  /// - Parameter masses: Mass of each particle, in yoctograms.
  /// - Parameter temperature: Temperature, in Kelvin.
  /// - Returns: Velocity of each particle, in nanometers per picosecond.
  ///
  /// Particles with zero mass should have zero velocity.
  static func boltzmannVelocities(
    masses: [Float],
    temperature: Float
  ) -> [SIMD3<Float>] {
    func createSystem() -> OpenMM_System {
      let system = OpenMM_System()
      for massInYg in masses {
        let massInAmu = massInYg * Float(MM4AmuPerYg)
        system.addParticle(mass: Double(massInAmu))
      }
      return system
    }
    let system = createSystem()
    
    func createEmptyPositions() -> OpenMM_Vec3Array {
      let output = OpenMM_Vec3Array(size: 0)
      for _ in masses.indices {
        let position: SIMD3<Double> = .zero
        output.append(position)
      }
      return output
    }
    
    func createContext(system: OpenMM_System) -> OpenMM_Context {
      let integrator = OpenMM_VerletIntegrator(
        stepSize: 0)
      let context = OpenMM_Context(
        system: system,
        integrator: integrator)
      context.positions = createEmptyPositions()
      context.setVelocitiesToTemperature(Double(temperature))
      return context
    }
    let context = createContext(system: system)
    
    func createVelocities(state: OpenMM_State) -> OpenMM_Vec3Array {
      let velocities = state.velocities
      guard velocities.size == masses.count else {
        fatalError("OpenMM array had incorrect size.")
      }
      return velocities
    }
    let state = context.state(types: .velocities)
    let velocities = createVelocities(state: state)
    
    var output: [SIMD3<Float>] = []
    for particleID in masses.indices {
      let velocity64 = velocities[particleID]
      let velocity32 = SIMD3<Float>(velocity64)
      output.append(velocity32)
    }
    
    // Strange crash when the C++ object for the OpenMM::State is deleted.
    withExtendedLifetime(state) { }
    return output
  }
  
  /// Perform velocity rescaling to exactly match the desired temperature.
  /// - Parameter velocities: Velocity of each particle, in nanometers per
  ///   picosecond.
  /// - Parameter masses: Mass of each particle, in yoctograms.
  /// - Parameter temperature: Temperature, in Kelvin.
  static func rescale(
    velocities: inout [SIMD3<Float>],
    masses: [Float],
    temperature: Float
  ) {
    var actualTotalKinetic: Double = .zero
    var expectedTotalKinetic: Double = .zero
    for particleID in masses.indices {
      let mass = masses[particleID]
      let velocity = velocities[particleID]
      let kineticEnergy = 0.5 * mass * (velocity * velocity).sum()
      
      let boltzmannInJPerK: Float = 1.380649e-23
      let boltzmannInZJPerK = boltzmannInJPerK * 1e21
      let thermalKineticEnergy = 1.5 * boltzmannInZJPerK * temperature
      
      if mass > 0 {
        actualTotalKinetic += Double(kineticEnergy)
        expectedTotalKinetic += Double(thermalKineticEnergy)
      }
    }
    
    var correctionFactor = Float(expectedTotalKinetic / actualTotalKinetic)
    correctionFactor = correctionFactor.squareRoot()
    
    for particleID in masses.indices {
      var velocity = velocities[particleID]
      velocity *= correctionFactor
      velocities[particleID] = velocity
    }
  }
  
  static func temperature(
    velocities: [SIMD3<Float>],
    masses: [Float]
  ) -> Double {
    var totalKineticEnergy: Double = .zero
    var degreesOfFreedom: Int = .zero
    
    for particleID in masses.indices {
      let mass = masses[particleID]
      let velocity = velocities[particleID]
      let kineticEnergy = 0.5 * mass * (velocity * velocity).sum()
      
      if mass > 0 {
        totalKineticEnergy += Double(kineticEnergy)
        degreesOfFreedom += 3
      }
    }
    
    // E = 1.5 * n * k * T
    // E = 0.5 * DOFs * k * T
    // T = E / (0.5 * DOFs * k)
    let boltzmannInJPerK: Double = 1.380649e-23
    let boltzmannInZJPerK = boltzmannInJPerK * 1e21
    var temperature = totalKineticEnergy
    temperature /= 0.5 * Double(degreesOfFreedom) * boltzmannInZJPerK
    return temperature
  }
}
