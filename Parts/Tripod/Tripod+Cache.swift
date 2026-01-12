import Foundation
import HDL

extension Tripod {
  // Procedure for generating a unique identifier for the current state.
  static func createKey(atoms: [Atom]) -> String {
    let key = Serialization.hash(atoms: atoms)
    
    // RFC 3548 encoding: https://www.rfc-editor.org/rfc/rfc3548#page-6
    // "/" -> "_"
    // "+" -> "-"
    var base64Key = key.base64EncodedString()
    do {
      // Fetch the null-terminated C string.
      var cString = base64Key.utf8CString
      for characterID in cString.indices {
        let byte = cString[characterID]
        let scalar = UnicodeScalar(UInt32(byte))!
        var character = Character(scalar)
        
        if character == "/" {
          character = "_"
        } else if character == "+" {
          character = "-"
        }
        cString[characterID] = CChar(character.asciiValue!)
      }
      base64Key = cString.withUnsafeBufferPointer {
        return String(cString: $0.baseAddress!)
      }
    }
    return base64Key
  }
  
  static func loadCachedTrajectory(
    atoms: [Atom]
  ) -> [[Atom]] {
    // Find the path.
    let packagePath = FileManager.default.currentDirectoryPath
    let key = createKey(atoms: atoms)
    
    // Load the cached trajectory.
    func loadData() -> Data {
      let fileURL = URL(
        filePath: "\(packagePath)/.build/Tripod/\(key)")
      if let data = try? Data(contentsOf: fileURL) {
        return data
      }
      
      let directoryURL = URL(
        filePath: "\(packagePath)/.build/Tripod")
      try! FileManager.default.createDirectory(
        at: directoryURL,
        withIntermediateDirectories: true)
      
      let originalFrames = runMinimization(atoms: atoms)
      let data = Serialization.encode(frames: originalFrames)
      
      let succeeded = FileManager.default.createFile(
        atPath: "\(packagePath)/.build/Tripod/\(key)",
        contents: data)
      guard succeeded else {
        fatalError("Could not create file.")
      }
      
      return data
    }
    
    // Ensure both branches return the exact same data. No discrepancies
    // on different program runs.
    let data = loadData()
    return Serialization.decode(frames: data)
  }
}

