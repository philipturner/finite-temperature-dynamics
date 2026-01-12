import Foundation
import GIFModule
import HDL
import MM4
import MolecularRenderer
import QuaternionModule
import xTB

// MARK: - User-Facing Options

let renderingOffline: Bool = false

// The simulation time per frame, in picoseconds. Frames are recorded and
// nominally played back at 60 FPS.
let frameSimulationTime: Double = 0.25 / 60
let simulationFrameCount: Int = 60 * 2 // 5 seconds for final render
let gifFrameSkipRate: Int = 1

// separate top-level constant for gifFrameCount
let freezeTimestamp: Float = 1.0
let gifFrameCount = 60 + simulationFrameCount + 0

// WARNING: Change the label between every video.
let tripodType: TripodType = .azatrane(.tin)
let feedstockType: FeedstockType = .carbonDimer(.bromine)
let videoLabel = "N3Sn-CCBr"

let temperature: Float = 300

// MARK: - Setup

var tripodDesc = TripodDescriptor()
tripodDesc.tripodType = tripodType
tripodDesc.feedstockType = feedstockType
var tripod = Tripod(descriptor: tripodDesc)
tripod.minimize()
let minimizedAtoms = tripod.cage.atoms + tripod.feedstock

let masses = minimizedAtoms
  .map(\.atomicNumber)
  .map(VerletIntegrator.mass(atomicNumber:))
var velocities = VerletIntegrator.boltzmannVelocities(
  masses: masses,
  temperature: temperature * 2)
VerletIntegrator.rescale(
  velocities: &velocities,
  masses: masses,
  temperature: temperature * 2)

var integratorDesc = VerletIntegratorDescriptor()
integratorDesc.masses = masses
var integrator = VerletIntegrator(descriptor: integratorDesc)
integrator.positions = minimizedAtoms.map(\.position)
integrator.velocities = velocities

xTB_Environment.verbosity = .muted
var calculatorDesc = xTB_CalculatorDescriptor()
calculatorDesc.atomicNumbers = minimizedAtoms.map(\.atomicNumber)
let calculator = xTB_Calculator(descriptor: calculatorDesc)

// MARK: - Simulation

@MainActor
func evaluateSinglepoint(positions: [SIMD3<Float>]) -> (
  forces: [SIMD3<Float>],
  energy: Double
) {
  calculator.molecule.positions = positions
  
  let forces = calculator.molecule.forces
  let energy = calculator.energy
  return (forces, energy)
}

@MainActor
func equilibriate() {
  let intervalDuration: Double = 0.010
  integrator.timeForConsole = nil
  
  for _ in 0..<10 {
    integrator.simulate(
      time: intervalDuration,
      singlepoint: evaluateSinglepoint)
    VerletIntegrator.rescale(
      velocities: &integrator.velocities,
      masses: masses,
      temperature: temperature)
  }
}
equilibriate()

var frames: [[Atom]] = []
@MainActor
func createFrame() -> [Atom] {
  var output: [Atom] = []
  for atomID in masses.indices {
    let atomicNumber = minimizedAtoms[atomID].atomicNumber
    let position = integrator.positions[atomID]
    let atom = Atom(position: position, atomicNumber: atomicNumber)
    output.append(atom)
  }
  return output
}
frames.append(createFrame())

for frameID in 1...simulationFrameCount {
  let time = Double(frameID) * frameSimulationTime
  integrator.timeForConsole = time
  integrator.simulate(
    time: frameSimulationTime,
    singlepoint: evaluateSinglepoint(positions:))
  frames.append(createFrame())
}

// MARK: - Launch Application

// Input: time in seconds
// Output: atoms
func interpolate(
  frames: [[Atom]],
  time: Float
) -> [Atom] {
  guard frames.count >= 1 else {
    fatalError("Need at least one frame to know size of atom list.")
  }
  
  let multiple60Hz = time * 60
  var lowFrame = Int(multiple60Hz.rounded(.down))
  var highFrame = lowFrame + 1
  var lowInterpolationFactor = Float(highFrame) - multiple60Hz
  var highInterpolationFactor = multiple60Hz - Float(lowFrame)
  
  if lowFrame < -1 {
    fatalError("This should never happen.")
  }
  if lowFrame >= frames.count - 1 {
    lowFrame = frames.count - 1
    highFrame = frames.count - 1
    lowInterpolationFactor = 1
    highInterpolationFactor = 0
  }
  
  var output: [Atom] = []
  for atomID in frames[0].indices {
    let lowAtom = frames[lowFrame][atomID]
    let highAtom = frames[highFrame][atomID]
    
    var position: SIMD3<Float> = .zero
    position += lowAtom.position * lowInterpolationFactor
    position += highAtom.position * highInterpolationFactor
    
    var atom = lowAtom
    atom.position = position
    output.append(atom)
  }
  return output
}

@MainActor
func createApplication() -> Application {
  // Set up the device.
  var deviceDesc = DeviceDescriptor()
  deviceDesc.deviceID = Device.fastestDeviceID
  let device = Device(descriptor: deviceDesc)
  
  // Set up the display.
  var displayDesc = DisplayDescriptor()
  displayDesc.device = device
  if renderingOffline {
    displayDesc.frameBufferSize = SIMD2<Int>(720, 1280)
  } else {
    displayDesc.frameBufferSize = SIMD2<Int>(1080, 1920)
  }
  if !renderingOffline {
    displayDesc.monitorID = device.fastestMonitorID
  }
  let display = Display(descriptor: displayDesc)
  
  // Set up the application.
  var applicationDesc = ApplicationDescriptor()
  applicationDesc.device = device
  applicationDesc.display = display
  if renderingOffline {
    applicationDesc.upscaleFactor = 1
  } else {
    applicationDesc.upscaleFactor = 3
  }
  
  applicationDesc.addressSpaceSize = 4_000_000
  applicationDesc.voxelAllocationSize = 500_000_000
  applicationDesc.worldDimension = 64
  let application = Application(descriptor: applicationDesc)
  
  return application
}
let application = createApplication()

@MainActor
func createTime() -> Float {
  if renderingOffline {
    let elapsedFrames = gifFrameSkipRate * application.frameID
    let frameRate: Int = 60
    let seconds = Float(elapsedFrames) / Float(frameRate)
    return seconds
  } else {
    let elapsedFrames = application.clock.frames
    let frameRate = application.display.frameRate
    let seconds = Float(elapsedFrames) / Float(frameRate)
    return seconds
  }
}

@MainActor
func modifyAtoms() {
  // Add one second of delay before and after rendering atoms, to just
  // render rotations. If not needed, we can cut it out during post-processing
  // with DaVinci Resolve.
  //
  // No special branches needed here, just add 60 frames to the GIF frame
  // count.
  let time = createTime()
  
  if time < freezeTimestamp {
    let atoms = frames.first!
    for atomID in atoms.indices {
      let atom = atoms[atomID]
      application.atoms[atomID] = atom
    }
  } else {
    let atoms = interpolate(
      frames: frames,
      time: time - freezeTimestamp)
    for atomID in atoms.indices {
      let atom = atoms[atomID]
      application.atoms[atomID] = atom
    }
  }
}

@MainActor
func modifyCamera() {
  let time = createTime()
  let angleDegrees = 0.1 * time * 360
  let timeRotation = Quaternion<Float>(
    angle: Float.pi / 180 * angleDegrees,
    axis: SIMD3(0, 1, 0))
  
  // Rotate the camera over time
  let focalPoint = SIMD3<Float>(0, 0, 0)
  var rotation = Quaternion<Float>(
    angle: Float.pi / 180 * 20,
    axis: SIMD3(-1, 0, 0))
  rotation = timeRotation * rotation
  let cameraDistance: Float = 5
  
  func rotate(_ vector: SIMD3<Float>) -> SIMD3<Float> {
    var output = rotation.act(on: vector)
    
    // Fix source of rendering error. Now the limiting factor is probably
    // internal to the renderer itself. Unsure exactly what's happening at
    // large distances.
    output /= (output * output).sum().squareRoot()
    return output
  }
  application.camera.basis.0 = rotate(SIMD3(1, 0, 0))
  application.camera.basis.1 = rotate(SIMD3(0, 1, 0))
  application.camera.basis.2 = rotate(SIMD3(0, 0, 1))
  application.camera.fovAngleVertical = Float.pi / 180 * 30
  
  var position = focalPoint
  position += rotation.act(on: SIMD3(0, 0, cameraDistance))
  application.camera.position = position
}

// Enter the run loop.
if !renderingOffline {
  application.run {
    modifyAtoms()
    modifyCamera()
    
    var image = application.render()
    image = application.upscale(image: image)
    application.present(image: image)
  }
} else {
  let frameBufferSize = application.display.frameBufferSize
  var gif = GIF(
    width: frameBufferSize[0],
    height: frameBufferSize[1])
  
  // Overall latency summary for offline mode:
  //
  // throughput @ 1440x1080, 60 FPS
  // macOS: 22.8 minutes / minute of content
  // Windows: 31.3 minutes / minute of content
  //
  // throughput @ 1280x720, 60 FPS
  // macOS: 13.5 minutes / minute of content
  // Windows: 18.5 minutes / minute of content
  //
  // Costs are probably agnostic to level of detail in the scene. On macOS, the
  // encoding latency was identical for an accidentally 100% black image.
  print("rendering frames")
  for _ in 0..<(gifFrameCount / gifFrameSkipRate) {
    let loopStartCheckpoint = Date()
    modifyAtoms()
    modifyCamera()
    
    // GPU-side bottleneck
    // throughput @ 1440x1080, 64 AO samples
    // macOS: 14-18 ms/frame
    // Windows: 50-70 ms/frame
    let image = application.render()
    
    // single-threaded bottleneck
    // throughput @ 1440x1080
    // macOS: 5 ms/frame
    // Windows: 47 ms/frame
    var gifImage = GIFModule.Image(
      width: frameBufferSize[0],
      height: frameBufferSize[1])
    for y in 0..<frameBufferSize[1] {
      for x in 0..<frameBufferSize[0] {
        let address = y * frameBufferSize[0] + x
        
        // Leaving this in the original SIMD4<Float16> causes a CPU-side
        // bottleneck on Windows.
        let pixel = SIMD4<Float>(image.pixels[address])
        
        // Don't clamp to [0, 255] range to avoid a minor CPU-side bottleneck.
        // It theoretically should never go outside this range; we just lose
        // the ability to assert this.
        let scaled = pixel * 255
        
        // On the Windows machine, '.toNearestOrEven' causes a massive
        // CPU-side bottleneck.
        let rounded = (scaled + 0.5).rounded(.down)
        
        // Avoid massive CPU-side bottleneck for unknown reason when casting
        // floating point vector to integer vector.
        let r = UInt8(rounded[0])
        let g = UInt8(rounded[1])
        let b = UInt8(rounded[2])
        
        let color = Color(
          red: r,
          green: g,
          blue: b)
        
        gifImage[y, x] = color
      }
    }
    
    // single-threaded bottleneck
    // throughput @ 1440x1080
    // macOS: 76 ms/frame
    // Windows: 271 ms/frame
    let quantization = OctreeQuantization(fromImage: gifImage)
    
    // For some reason, DaVinci Resolve imports 20 FPS clips as 25 FPS. So I
    // change delayTime to 4 when exporting to DaVinci Resolve.
    let frame = Frame(
      image: gifImage,
      delayTime: 4,
      localQuantization: quantization)
    gif.frames.append(frame)
    
    let loopEndCheckpoint = Date()
    print(loopEndCheckpoint.timeIntervalSince(loopStartCheckpoint))
  }
  
  // multi-threaded bottleneck
  // throughput @ 1440x1080
  // macOS: 252 ms/frame
  // Windows: 174 ms/frame (abnormally fast compared to macOS)
  print("encoding GIF")
  let encodeStartCheckpoint = Date()
  let data = try! gif.encoded()
  let encodeEndCheckpoint = Date()
  
  let encodedSizeRepr = String(format: "%.1f", Float(data.count) / 1e6)
  print("encoded size:", encodedSizeRepr, "MB")
  print(encodeEndCheckpoint.timeIntervalSince(encodeStartCheckpoint))
  
  // SSD access bottleneck
  //
  // latency @ 1440x1080, 10 frames, 2.1 MB
  // macOS: 1.6 ms
  // Windows: 16.3 ms
  //
  // latency @ 1440x1080, 60 frames, 12.4 MB
  // macOS: 4.1 ms
  // Windows: 57.7 ms
  //
  // Order of magnitude, 1 minute of video is 1 GB of GIF.
  let packagePath = FileManager.default.currentDirectoryPath
  let filePath = "\(packagePath)/.build/\(videoLabel).gif"
  let succeeded = FileManager.default.createFile(
    atPath: filePath,
    contents: data)
  guard succeeded else {
    fatalError("Could not write to file.")
  }
}
