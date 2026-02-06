use chemfiles::Trajectory;
use molly::XTCReader;

fn round_to(v: f32, decimals: u32) -> f32 {
    let d = f32::powi(10.0, decimals as i32);
    (v * d).round() / d
}

fn main() -> std::io::Result<()> {
    let mut args = std::env::args().skip(1);
    let path = args.next().expect("please provide one xtc trajectory path");
    let decimals: u32 = args.next().and_then(|d| d.parse().ok()).unwrap_or(3);
    dbg!(decimals);

    let file = std::fs::File::open(&path)?;
    let mut reader = XTCReader::new(file);
    let mut trajectory = Trajectory::open(&path, 'r').unwrap();
    let mut cfframe = chemfiles::Frame::new();
    let mut n = 0;
    let mut natoms = 0;
    let mut frame = molly::Frame::default();
    while reader.read_frame(&mut frame).is_ok() {
        trajectory.read(&mut cfframe).unwrap();

        for (a, &b) in frame.coords().zip(cfframe.positions()) {
            let a = a.map(|v| round_to(v, decimals));
            let b = b.map(|v| v as f32 * 0.1).map(|v| round_to(v, decimals));
            // eprintln!("a, b = {a:?}\t\t{b:?}");
            assert_eq!(a, b);
        }

        natoms = frame.positions.len() / 3;
        n += 1;
    }
    eprintln!("compare: read {n} frames");
    assert_eq!(natoms, cfframe.positions().len());
    eprintln!("{} atoms", cfframe.positions().len());

    Ok(())
}
