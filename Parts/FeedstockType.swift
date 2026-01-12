import HDL
import QuaternionModule

enum FeedstockType {
  case hydrogen
  case carbon(Element?, Element?, Element?)
  case silicon(Element?, Element?, Element?)
  case germanium(Element?, Element?, Element?)
  case carbonDimer(Element?)
  
  // Count the elements in the optional parameters, create a compacted list
  // without changing the original order.
  //
  // carbon:
  // case 1 (carbene)     120° zenith, no azimuth
  // case 2 (methylene)   120° zenith, 180° azimuth
  // case 3 (methane)     109.5° zenith, 120° azimuth
  //
  // silicon, germanium:
  // case 1 (silene)      109.5° zenith, no azimuth
  // case 2 (silylene)    109.5° zenith, 120° azimuth
  // case 2 (silane)      109.5° zenith, 120° azimuth
  //
  // base atom is placed at (0, 0, 0) and first in the list
  var atoms: [Atom] {
    func createPassivatorList(
      _ passivator1: Element?,
      _ passivator2: Element?,
      _ passivator3: Element?
    ) -> [Element] {
      var output: [Element] = []
      if let passivator1 {
        output.append(passivator1)
      }
      if let passivator2 {
        output.append(passivator2)
      }
      if let passivator3 {
        output.append(passivator3)
      }
      return output
    }
    
    switch self {
    case .hydrogen:
      return [Atom(position: .zero, element: .hydrogen)]
      
    case .carbon(let passivator1, let passivator2, let passivator3):
      let passivators = createPassivatorList(
        passivator1,
        passivator2,
        passivator3)
      if passivators.count <= 2 {
        return Self.createMethyleneAtoms(
          center: .carbon,
          passivators: passivators)
      } else {
        return Self.createMethaneAtoms(
          center: .carbon,
          passivators: passivators)
      }
      
    case .silicon(let passivator1, let passivator2, let passivator3):
      let passivators = createPassivatorList(
        passivator1,
        passivator2,
        passivator3)
      return Self.createMethaneAtoms(
        center: .silicon,
        passivators: passivators)
      
    case .germanium(let passivator1, let passivator2, let passivator3):
      let passivators = createPassivatorList(
        passivator1,
        passivator2,
        passivator3)
      return Self.createMethaneAtoms(
        center: .germanium,
        passivators: passivators)
      
    case .carbonDimer(let passivator):
      return Self.createCarbonDimerAtoms(
        passivator: passivator)
    }
  }
  
  // Abbreviated label for the debug string representation.
  //
  // Examples: H, CHH, CHHBr, SiHHH, CC, CCH
  var description: String {
    func repr(passivator: Element) -> String {
      switch passivator {
      case .hydrogen:
        return "H"
      case .fluorine:
        return "F"
      case .chlorine:
        return "Cl"
      case .bromine:
        return "Br"
      default:
        fatalError("Unsupported passivator element: \(passivator)")
      }
    }
    
    func repr(
      _ passivator1: Element?,
      _ passivator2: Element?,
      _ passivator3: Element?
    ) -> String {
      var output: String = ""
      if let passivator1 {
        output += repr(passivator: passivator1)
      }
      if let passivator2 {
        output += repr(passivator: passivator2)
      }
      if let passivator3 {
        output += repr(passivator: passivator3)
      }
      return output
    }
    
    switch self {
    case .hydrogen:
      return "H"
    case .carbon(let passivator1, let passivator2, let passivator3):
      var output = "C"
      output += repr(passivator1, passivator2, passivator3)
      return output
    case .silicon(let passivator1, let passivator2, let passivator3):
      var output = "Si"
      output += repr(passivator1, passivator2, passivator3)
      return output
    case .germanium(let passivator1, let passivator2, let passivator3):
      var output = "Ge"
      output += repr(passivator1, passivator2, passivator3)
      return output
    case .carbonDimer(let passivator):
      var output = "CC"
      if let passivator {
        output += repr(passivator: passivator)
      }
      return output
    }
  }
}

extension FeedstockType {
  static func createMethyleneAtoms(
    center centerElement: Element,
    passivators passivatorElements: [Element]
  ) -> [Atom] {
    guard passivatorElements.count >= 1,
          passivatorElements.count <= 2 else {
      fatalError("Unexpected number of passivators.")
    }
    
    let center = Atom(
      position: SIMD3.zero,
      element: centerElement)
    let zenithRotation = Quaternion<Float>(
      angle: Float.pi / 180 * -120,
      axis: SIMD3(1, 0, 0))
    
    var passivators: [Atom] = []
    for passivatorID in passivatorElements.indices {
      let passivatorElement = passivatorElements[passivatorID]
      let bondLength =
      centerElement.covalentRadius + passivatorElement.covalentRadius
      
      var position = SIMD3<Float>(0, -bondLength, 0)
      position = zenithRotation.act(on: position)
      
      let azimuthRotation = Quaternion<Float>(
        angle: Float.pi / 180 * 180 * Float(passivatorID),
        axis: SIMD3(0, 1, 0))
      position = azimuthRotation.act(on: position)
      
      let passivator = Atom(
        position: position,
        element: passivatorElement)
      passivators.append(passivator)
    }
    return [center] + passivators
  }
  
  static func createMethaneAtoms(
    center centerElement: Element,
    passivators passivatorElements: [Element]
  ) -> [Atom] {
    guard passivatorElements.count >= 1,
          passivatorElements.count <= 3 else {
      fatalError("Unexpected number of passivators.")
    }
    
    let center = Atom(
      position: SIMD3.zero,
      element: centerElement)
    let zenithRotation = Quaternion<Float>(
      angle: Float.pi / 180 * -109.5,
      axis: SIMD3(1, 0, 0))
    
    var passivators: [Atom] = []
    for passivatorID in passivatorElements.indices {
      let passivatorElement = passivatorElements[passivatorID]
      let bondLength =
      centerElement.covalentRadius + passivatorElement.covalentRadius
      
      var position = SIMD3<Float>(0, -bondLength, 0)
      position = zenithRotation.act(on: position)
      
      let azimuthRotation = Quaternion<Float>(
        angle: Float.pi / 180 * 120 * Float(passivatorID),
        axis: SIMD3(0, 1, 0))
      position = azimuthRotation.act(on: position)
      
      let passivator = Atom(
        position: position,
        element: passivatorElement)
      passivators.append(passivator)
    }
    return [center] + passivators
  }
  
  static func createCarbonDimerAtoms(
    passivator passivatorElement: Element?
  ) -> [Atom] {
    let carbon1 = Atom(
      position: SIMD3.zero,
      element: .carbon)
    let carbon2 = Atom(
      position: SIMD3(0, 0.120, 0),
      element: .carbon)
    guard let passivatorElement else {
      return [carbon1, carbon2]
    }
    
    let bondLength =
    Element.carbon.covalentRadius + passivatorElement.covalentRadius
    
    let position = carbon2.position + SIMD3(0, bondLength, 0)
    let passivator = Atom(
      position: position,
      element: passivatorElement)
    
    return [carbon1, carbon2, passivator]
  }
}
