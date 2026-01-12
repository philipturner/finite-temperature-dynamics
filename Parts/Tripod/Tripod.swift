import HDL
import xTB

struct TripodDescriptor {
  /// Required.
  var tripodType: TripodType?
  
  /// Optional.
  var feedstockType: FeedstockType?
}

struct Tripod {
  var cage: Topology
  var feedstock: [Atom]
  
  init(descriptor: TripodDescriptor) {
    guard let tripodType = descriptor.tripodType else {
      fatalError("Descriptor was incomplete.")
    }
    self.cage = tripodType.topology
    
    if let feedstockType = descriptor.feedstockType {
      self.feedstock = feedstockType.atoms
      
      // Shift the atoms up by the bond length.
      let apexAtom = cage.atoms[0]
      let apexAtomRadius = apexAtom.element.covalentRadius
      let centerAtomRadius = feedstock[0].element.covalentRadius
      let bondLength = apexAtomRadius + centerAtomRadius
      
      for atomID in feedstock.indices {
        var atom = feedstock[atomID]
        atom.position.y += bondLength
        feedstock[atomID] = atom
      }
    } else {
      self.feedstock = []
    }
  }
}
