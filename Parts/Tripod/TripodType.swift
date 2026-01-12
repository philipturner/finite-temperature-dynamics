import HDL
import QuaternionModule

enum TripodType {
  case adamantane(Element)
  case carbatrane(Element)
  case azatrane(Element)
  
  var topology: Topology {
    let skeleton: Topology
    
    switch self {
    case .adamantane(let apex):
      let allowedApexElements: [Element] = [
        .carbon, .silicon, .germanium
      ]
      guard allowedApexElements.contains(apex) else {
        fatalError("Invalid apex for cage type.")
      }
      skeleton = Self.createAdamantaneSkeleton(
        apex: apex)
      
    case .carbatrane(let apex):
      let allowedApexElements: [Element] = [
        .silicon, .germanium, .tin
      ]
      guard allowedApexElements.contains(apex) else {
        fatalError("Invalid apex for cage type.")
      }
      skeleton = Self.createAtraneSkeleton(
        apex: apex,
        isAzatrane: false)
      
    case .azatrane(let apex):
      let allowedApexElements: [Element] = [
        .silicon, .germanium, .tin
      ]
      guard allowedApexElements.contains(apex) else {
        fatalError("Invalid apex for cage type.")
      }
      skeleton = Self.createAtraneSkeleton(
        apex: apex,
        isAzatrane: true)
    }
    
    func createHydrogen(
      atomID: UInt32,
      orbital: SIMD3<Float>
    ) -> Atom {
      let atom = skeleton.atoms[Int(atomID)]
      
      var bondLength = atom.element.covalentRadius
      bondLength += Element.hydrogen.covalentRadius
      
      let position = atom.position + bondLength * orbital
      return Atom(position: position, element: .hydrogen)
    }
    
    let orbitalLists = skeleton.nonbondingOrbitals()
    
    var insertedAtoms: [Atom] = []
    var insertedBonds: [SIMD2<UInt32>] = []
    for atomID in skeleton.atoms.indices {
      let atom = skeleton.atoms[atomID]
      guard atom.atomicNumber == 6 else {
        continue
      }
      
      // Skip the apex atom.
      if atomID == 0 {
        continue
      }
      
      let orbitalList = orbitalLists[atomID]
      for orbital in orbitalList {
        let hydrogen = createHydrogen(
          atomID: UInt32(atomID),
          orbital: orbital)
        let hydrogenID = skeleton.atoms.count + insertedAtoms.count
        insertedAtoms.append(hydrogen)
        
        let bond = SIMD2(
          UInt32(atomID),
          UInt32(hydrogenID))
        insertedBonds.append(bond)
      }
    }
    
    var output = skeleton
    output.atoms += insertedAtoms
    output.bonds += insertedBonds
    return output
  }
}

extension TripodType {
  // The way the lattice compiles, the apex atom ends up first in the list.
  static func createAdamantaneSkeleton(
    apex apexElement: Element,
  ) -> Topology {
    let lattice = Lattice<Cubic> { h, k, l in
      Bounds { 1 * (h + k + l) }
      Material { .elemental(.carbon) }
      
      Volume {
        Origin { 0.5 * (h + k + l) }
        Concave {
          Plane { -h }
          Plane { -k }
          Plane { -l }
        }
        Replace { .atom(apexElement) }
      }
    }
    
    var reconstruction = Reconstruction()
    reconstruction.atoms = lattice.atoms
    reconstruction.material = .elemental(.carbon)
    var topology = reconstruction.compile()
    
    for atomID in topology.atoms.indices {
      var atom = topology.atoms[atomID]
      
      // Shift all atoms by -0.25 * latticeConstant
      let latticeConstant = Constant(.square) {
        .elemental(.carbon)
      }
      atom.position -= 0.25 * latticeConstant
      
      // Rotate the basis
      let basisX = SIMD3<Float>(-1, 0, 1) / Float(2).squareRoot()
      let basisY = SIMD3<Float>(-1, -1, -1) / Float(3).squareRoot()
      let basisZ = SIMD3<Float>(1, -2, 1) / Float(6).squareRoot()
      let newX = (atom.position * basisX).sum()
      let newY = (atom.position * basisY).sum()
      let newZ = (atom.position * basisZ).sum()
      atom.position = SIMD3(newX, newY, newZ)
      
      topology.atoms[atomID] = atom
    }
    
    // Insert the thiol linkers.
    do {
      let orbitalLists = topology.nonbondingOrbitals()
      
      var insertedAtoms: [Atom] = []
      var insertedBonds: [SIMD2<UInt32>] = []
      for atomID in topology.atoms.indices {
        let atom = topology.atoms[atomID]
        let orbitalList = orbitalLists[atomID]
        
        // Only select atoms in the bridgehead position.
        guard orbitalList.count == 1 else {
          continue
        }
        guard atom.position.y < -0.05 else {
          continue
        }
        
        func ccBondLength() -> Float {
          2 * Element.carbon.covalentRadius
        }
        func csBondLength() -> Float {
          Element.carbon.covalentRadius + Element.sulfur.covalentRadius
        }
        func shBondLength() -> Float {
          Element.sulfur.covalentRadius + Element.hydrogen.covalentRadius
        }
        func createThiolHDirection() -> SIMD3<Float> {
          var output = orbitalList.first!
          output.y = 0
          output /= (output * output).sum().squareRoot()
          return output
        }
        
        var carbon0Position = atom.position
        carbon0Position += orbitalList.first! * ccBondLength()
        let carbon0 = Atom(position: carbon0Position, element: .carbon)
        
        var sulfur1Position = carbon0.position
        sulfur1Position += SIMD3<Float>(0, -1, 0) * csBondLength()
        let sulfur1 = Atom(position: sulfur1Position, element: .sulfur)
        
        var hydrogen2Position = sulfur1.position
        hydrogen2Position += createThiolHDirection() * shBondLength()
        let hydrogen2 = Atom(position: hydrogen2Position, element: .hydrogen)
        
        // The base ID for inserted atoms to start at.
        let baseAtomID = topology.atoms.count + insertedAtoms.count
        insertedAtoms.append(carbon0)
        insertedAtoms.append(sulfur1)
        insertedAtoms.append(hydrogen2)
        
        insertedBonds.append(
          SIMD2(UInt32(atomID), UInt32(baseAtomID + 0)))
        insertedBonds.append(
          SIMD2(UInt32(baseAtomID + 0), UInt32(baseAtomID + 1)))
        insertedBonds.append(
          SIMD2(UInt32(baseAtomID + 1), UInt32(baseAtomID + 2)))
      }
      topology.atoms += insertedAtoms
      topology.bonds += insertedBonds
    }
    
    return topology
  }
  
  static func createAtraneSkeleton(
    apex apexElement: Element,
    isAzatrane: Bool
  ) -> Topology {
    var topology = Topology()
    topology.atoms += [
      Atom(position: SIMD3(0.00, 0.00, 0.00), element: apexElement),
      Atom(position: SIMD3(0.00, -0.26, -0.00), element: .nitrogen),
    ]
    
    for legID in 0..<3 {
      let baseAtomID = topology.atoms.count
      var insertedAtoms: [Atom] = []
      var insertedBonds: [SIMD2<UInt32>] = []
      
      // Isolate the temporary variables for this part in a contained scope.
      do {
        let carbon0 = Atom(
          position: SIMD3(0.15, -0.28, 0.00), element: .carbon)
        insertedAtoms.append(carbon0)
        insertedBonds.append(
          SIMD2(UInt32(1), UInt32(baseAtomID + 0)))
        
        let carbon1 = Atom(
          position: SIMD3(0.23, -0.15, -0.05), element: .carbon)
        insertedAtoms.append(carbon1)
        insertedBonds.append(
          SIMD2(UInt32(baseAtomID + 0), UInt32(baseAtomID + 1)))
        
        let nitrogen2 = Atom(
          position: SIMD3(0.23, -0.00, 0.00),
          element: isAzatrane ? .nitrogen : .carbon)
        insertedAtoms.append(nitrogen2)
        insertedBonds.append(
          SIMD2(UInt32(baseAtomID + 1), UInt32(baseAtomID + 2)))
        insertedBonds.append(
          SIMD2(UInt32(baseAtomID + 2), UInt32(0)))
        
        let carbon3 = Atom(
          position: SIMD3(0.38, -0.20, -0.05), element: .carbon)
        insertedAtoms.append(carbon3)
        insertedBonds.append(
          SIMD2(UInt32(baseAtomID + 1), UInt32(baseAtomID + 3)))
        
        let sulfur4 = Atom(
          position: SIMD3(0.45, -0.36, -0.05), element: .sulfur)
        insertedAtoms.append(sulfur4)
        insertedBonds.append(
          SIMD2(UInt32(baseAtomID + 3), UInt32(baseAtomID + 4)))
        
        let hydrogen5 = Atom(
          position: SIMD3(0.57, -0.36, -0.05), element: .hydrogen)
        insertedAtoms.append(hydrogen5)
        insertedBonds.append(
          SIMD2(UInt32(baseAtomID + 4), UInt32(baseAtomID + 5)))
      }
      
      if isAzatrane {
        let carbon6 = Atom(
          position: SIMD3(0.32, 0.08, -0.05), element: .carbon)
        insertedAtoms.append(carbon6)
        insertedBonds.append(
          SIMD2(UInt32(baseAtomID + 2), UInt32(baseAtomID + 6)))
        
        let hydrogen7 = Atom(
          position: SIMD3(0.32, 0.20, -0.05), element: .hydrogen)
        insertedAtoms.append(hydrogen7)
        insertedBonds.append(
          SIMD2(UInt32(baseAtomID + 6), UInt32(baseAtomID + 7)))
      }
      
      // Apply the rotation transform to all atoms, just before inserting.
      let angleDegrees = Float(legID) * 120
      let rotation = Quaternion<Float>(
        angle: Float.pi / 180 * angleDegrees,
        axis: SIMD3(0, 1, 0))
      for relativeAtomID in insertedAtoms.indices {
        var atom = insertedAtoms[relativeAtomID]
        atom.position = rotation.act(on: atom.position)
        insertedAtoms[relativeAtomID] = atom
      }
      topology.atoms += insertedAtoms
      topology.bonds += insertedBonds
    }
    
    return topology
  }
}
