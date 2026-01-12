import HDL
import var MM4.MM4YgPerAmu
import xTB

extension Tripod {
  // Replaces the atoms with the minimized structure.
  //
  // Recycles a previous minimization on disk (if availiable). To properly use
  // the cache, only call minimize() once after generating the compiled
  // structure.
  //
  // Returns all frames of the minimization.
  @discardableResult
  mutating func minimize() -> [[Atom]] {
    var compiledAtoms: [Atom] = []
    compiledAtoms += cage.atoms
    compiledAtoms += feedstock
    
    let trajectory = Self.loadCachedTrajectory(atoms: compiledAtoms)
    let minimizedAtoms = trajectory.last!
    
    do {
      let baseAddress: Int = 0
      for atomID in cage.atoms.indices {
        let atom = minimizedAtoms[baseAddress + atomID]
        cage.atoms[atomID] = atom
      }
    }
    
    do {
      let baseAddress: Int = cage.atoms.count
      for atomID in feedstock.indices {
        let atom = minimizedAtoms[baseAddress + atomID]
        feedstock[atomID] = atom
      }
    }
    
    return trajectory
  }
  
  static func runMinimization(atoms: [Atom]) -> [[Atom]] {
    xTB_Environment.verbosity = .muted
    defer { xTB_Environment.verbosity = .minimal }
    
    var calculatorDesc = xTB_CalculatorDescriptor()
    calculatorDesc.atomicNumbers = atoms.map(\.atomicNumber)
    let calculator = xTB_Calculator(descriptor: calculatorDesc)
    
    func createMinimization() -> FIREMinimization {
      var minimizationDesc = FIREMinimizationDescriptor()
      minimizationDesc.masses = atoms.map {
        if $0.atomicNumber == 1 {
          return Float(4.0 * MM4YgPerAmu)
        } else {
          return Float(12.011 * MM4YgPerAmu)
        }
      }
      minimizationDesc.positions = atoms.map(\.position)
      return FIREMinimization(descriptor: minimizationDesc)
    }
    var minimization = createMinimization()
    
    func createSulfurIDs() -> [UInt32] {
      var output: [UInt32] = []
      for atomID in atoms.indices {
        let atom = atoms[atomID]
        if atom.element == .sulfur {
          output.append(UInt32(atomID))
        }
      }
      guard output.count == 3 else {
        fatalError("Failed to locate all the sulfurs on the legs.")
      }
      return output
    }
    let sulfurIDs = createSulfurIDs()
    
    func createFrame() -> [Atom] {
      var output: [Atom] = []
      for atomID in atoms.indices {
        var atom = atoms[atomID]
        let position = minimization.positions[atomID]
        atom.position = position
        output.append(atom)
      }
      return output
    }
    var frames: [[Atom]] = []
    
    let maxFrameCount: Int = 1000
    print()
    for trialID in 0..<maxFrameCount {
      frames.append(createFrame())
      calculator.molecule.positions = minimization.positions
      
      // Enforce the constraints on leg sulfurs.
      var forces = calculator.molecule.forces
      do {
        var forceAccumulator: Float = .zero
        for atomID in sulfurIDs {
          let force = forces[Int(atomID)]
          forceAccumulator += force.y
        }
        forceAccumulator /= 3
        for atomID in sulfurIDs {
          var force = forces[Int(atomID)]
          force.y = forceAccumulator
          forces[Int(atomID)] = force
        }
      }
      
      var maximumForce: Float = .zero
      for atomID in calculator.molecule.atomicNumbers.indices {
        if minimization.anchors.contains(UInt32(atomID)) {
          continue
        }
        let force = forces[atomID]
        let forceMagnitude = (force * force).sum().squareRoot()
        maximumForce = max(maximumForce, forceMagnitude)
      }
      
      print("time: \(Format.time(minimization.time))", terminator: " | ")
      print("energy: \(Format.energy(calculator.energy))", terminator: " | ")
      print("max force: \(Format.force(maximumForce))", terminator: " | ")
      
      let converged = minimization.step(forces: forces)
      if !converged {
        print("Δt: \(Format.time(minimization.Δt))", terminator: " | ")
      }
      print()
      
      // Enforce the constraints on leg sulfurs.
      do {
        var positions = minimization.positions
        var positionAccumulator: Float = .zero
        for atomID in sulfurIDs {
          let position = positions[Int(atomID)]
          positionAccumulator += position.y
        }
        positionAccumulator /= 3
        for atomID in sulfurIDs {
          var position = positions[Int(atomID)]
          position.y = positionAccumulator
          positions[Int(atomID)] = position
        }
        minimization.positions = positions
      }
      
      if converged {
        frames.append(createFrame())
        break
      } else if trialID == maxFrameCount - 1 {
        // All structures should converge now. Change back to a silent warning
        // if things get really bad at some point.
        fatalError("WARNING: failed to converge!")
      }
    }
    
    return frames
  }
}
