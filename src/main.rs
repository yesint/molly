//! Filter an xtc trajectory, quickly.
//!
//! By Marieke Westendorp, 2024.
//! <ma3ke.cyber@gmail.com>
use std::fs::File;
use std::io::{self, BufWriter, Read, Seek, SeekFrom, Write};
use std::num::{NonZeroU64, ParseIntError};
use std::path::PathBuf;
use std::str::FromStr;

use clap::Parser;
use molly::buffer::{Buffer, UnBuffered};
use molly::reader::{read_nbytes, NBYTES_POSITIONS_PRELUDE};
use molly::selection::{AtomSelection, FrameSelection, Range};
use molly::{padding, read_positions, Frame, Header, Magic, XTCReader, XTC_1995_MAX_NATOMS};

fn filter_frames(
    reader: &mut XTCReader<File>,
    writer: &mut BufWriter<File>,
    args: WriteArgs,
) -> std::io::Result<()> {
    // Forced magic drinking.
    let forced_magic = args
        .force_magic
        .map(Magic::try_from)
        .transpose()
        .map_err(std::io::Error::other)?;

    let mut scratch = Vec::new();

    let frame_selection = args.frame_selection.unwrap_or_default();
    let atom_selection = args.atom_selection.unwrap_or_default();

    let until = if args.reverse || args.reverse_frame_selection {
        None
    } else {
        frame_selection.until()
    };
    let mut offsets = reader.determine_offsets(until)?;
    let mut range: Box<[usize]> = (0..offsets.len()).collect();
    let enumerated_offsets = {
        // Reversing the frame order and reversing the frame selection have some non-obvious
        // interplays.
        match (args.reverse, args.reverse_frame_selection) {
            (true, true) => {
                offsets.reverse();
            }
            (true, false) => {
                offsets.reverse();
                range.reverse();
            }
            (false, true) => {
                range.reverse();
            }
            (false, false) => {}
        }
        range.iter().zip(offsets.iter().copied())
    };

    let mut stdout = std::io::stdout();
    let mut frame = Frame::default();
    for (&idx, offset) in enumerated_offsets {
        match frame_selection.is_included(idx) {
            Some(true) => {}
            Some(false) => continue,
            // If we are reversed in some way, we can't just stop early.
            None if args.reverse || args.reverse_frame_selection => continue,
            None => break,
        }

        // Go to the start of this frame.
        reader.file.seek(SeekFrom::Start(offset))?;

        // Start of by reading the header.
        let header = reader.read_header()?;

        if args.times || args.steps {
            if args.times {
                write!(stdout, "{:.3}\t", header.time)?;
            }
            if args.steps {
                write!(stdout, "{}", header.step)?;
            }
            writeln!(stdout, "")?;
        }

        // Now, we read the atoms.
        let natoms_frame = header.natoms; // The number of atoms specified for the frame.
        let nbytes = if natoms_frame <= 9 {
            // In this case, the positions are uncompressed. Each consists of three f32s, so we're
            // done pretty quickly.
            reader.read_smol_positions(natoms_frame, &mut frame, &atom_selection)?
        } else {
            let nbytes = match args.is_buffered {
                false => read_positions::<UnBuffered, File>(
                    &mut reader.file,
                    natoms_frame,
                    &mut scratch,
                    &mut frame,
                    &atom_selection,
                    header.magic,
                )?,
                true => read_positions::<Buffer, File>(
                    &mut reader.file,
                    natoms_frame,
                    &mut scratch,
                    &mut frame,
                    &atom_selection,
                    header.magic,
                )?,
            };
            reader.step += 1;
            nbytes
        };

        // The number of atoms we are actually interested in for our output. Important to know
        // since it may be the case that more atoms are selected than are in the frame.
        let natoms = frame.natoms();
        // Reset to the start of the frame again, and skip the header.
        let offset_and_header = offset + Header::SIZE as u64;
        reader.file.seek(SeekFrom::Start(offset_and_header))?;

        // Redefine the header to reflect our changes.
        let old_magic = header.magic;
        let header = Header {
            magic: forced_magic.unwrap_or(header.magic),
            natoms,
            natoms_repeated: natoms,
            ..header
        };
        // And write it.
        writer.write_all(&header.to_be_bytes())?;

        if header.magic == Magic::Xtc1995 && header.natoms > XTC_1995_MAX_NATOMS {
            eprintln!(
                "WARNING: The number of atoms to be written out ({}) \
                    exceeds the maximum number of atoms for the {} magic number \
                    ({XTC_1995_MAX_NATOMS} atoms)",
                header.natoms,
                Magic::Xtc1995
            )
        }

        if natoms <= 9 {
            // The number of positions is small. We encode the positions as uncompressed floats.
            for pos in &frame.positions {
                writer.write_all(&pos.to_be_bytes())?;
            }
        } else {
            // TODO: Consider 're-using' the scratch buffer!! It will contain (more than) the bytes we want to write out!
            // TODO: Invent some sort of SCRATCH mechanism here again.

            // Just copy over the precision, prelude, followed by the section of compressed bytes.
            let mut precision = [0; 4];
            reader.file.read_exact(&mut precision)?;
            writer.write_all(&precision)?;

            // Copy over the prelude, since that remains exactly the same.
            let mut prelude = [0; NBYTES_POSITIONS_PRELUDE];
            reader.file.read_exact(&mut prelude)?;
            writer.write_all(&prelude)?;

            // Note that we need to read according to the `old_magic`, since that describes the
            // data that we are about to read from. This matters since the magic may have been
            // forced to a different value through the secret command line option :)
            let nbytes_old = read_nbytes(&mut reader.file, old_magic)?;
            // Check whether we totally messed up.
            assert!(
                nbytes <= nbytes_old as usize,
                "the new number of bytes ({nbytes}) must never be greater than the old number of bytes ({nbytes_old})"
            );

            // Write the new number of upcoming bytes.
            match header.magic {
                Magic::Xtc1995 => writer.write_all(&(nbytes as u32).to_be_bytes())?,
                Magic::Xtc2023 => writer.write_all(&(nbytes as u64).to_be_bytes())?,
            }
            // Note that we are dealing with xdr padding, here! (32-bit blocks.)
            let mut bytes = vec![0; nbytes + padding(nbytes)];
            reader.file.read_exact(&mut bytes[..nbytes])?;
            writer.write_all(&bytes)?;
        }
    }

    Ok(())
}

fn frame_selection_parser(selection: &str) -> Result<FrameSelection, ParseIntError> {
    let mut components = selection.split(':');
    let start = components
        .next()
        .and_then(|s| if s.is_empty() { None } else { Some(s) })
        .map(|s| s.parse())
        .transpose()?;
    let end = components
        .next()
        .and_then(|s| if s.is_empty() { None } else { Some(s) })
        .map(|s| s.parse())
        .transpose()?;
    let step = components
        .next()
        .and_then(|s| if s.is_empty() { None } else { Some(s) })
        .map(|s| NonZeroU64::from_str(s))
        .transpose()?;
    Ok(FrameSelection::Range(Range::new(start, end, step)))
}

fn atom_selection_parser(selection: &str) -> Result<AtomSelection, ParseIntError> {
    let until: u32 = selection.parse()?;
    Ok(AtomSelection::Until(until))
}

// TODO: Consider making this one of several subcommands. This one could be called something like
// `molly filter ...`. Another would be `molly info` or `molly summary` or something.
/// Filter an xtc trajectory according to frame and atom selections.
///
/// By Marieke Westendorp, 2024.
/// <ma3ke.cyber@gmail.com>
#[derive(Parser)]
#[command(version)]
struct Args {
    /// Input path (xtc).
    input: PathBuf,

    #[command(flatten)]
    write: Option<WriteArgs>,

    /// Print a summary of the trajectory to standard output and exit.
    ///
    /// Conflicts with any options for writing.
    ///
    /// Currently, selections have no effect on the info displayed.
    #[arg(long, conflicts_with = "WriteArgs")]
    info: bool,
}

#[derive(Parser)]
struct WriteArgs {
    /// Output path (xtc).
    output: PathBuf,

    /// Frame selection in the format `start:stop:step`. Each of these values optional.
    ///
    // TODO: Make these examples into unit tests for the frame_selection_parser and its atom counterpart.
    // TODO: Verify that I didn't make any mistakes in these examples, once everything is up and running.
    /// - `:100` will select the first 100 frames.
    ///
    /// - `3:14` will select the 4th up to and including the 14th frames, 11 frames in total.
    ///
    /// - `:100:2` will select every second frame from the the first 100 frames, 50 in total.
    #[arg(short, long, value_parser=frame_selection_parser)]
    frame_selection: Option<FrameSelection>,

    // TODO: Consider explaining why this seemingly silly limitation exists. It may be confusing
    // to just drop it here, but explaining it is also quite the ride.
    /// Atom selection single `stop` value.
    ///
    /// For each frame that is read, the compressed positions up to the provided index will be
    /// stored into the output file.
    ///
    // TODO: Verify that I didn't make any mistakes in these examples, once everything is up and running.
    /// - `1312` selects the first 1312 atoms.
    ///
    /// Note that according to the xtc format, when the number of atoms in the frame is equal to
    /// or less than 9 (natoms <= 9), the atoms will be stored in an uncompressed manner.
    #[arg(short, long, value_parser=atom_selection_parser)]
    atom_selection: Option<AtomSelection>,

    // TODO: Add some {on, off, auto} enum?
    /// Use non-buffered reading mode. (Reading mode is buffered by default.)
    ///
    /// This may be faster under some circumstances, especially when the atom selection includes
    /// most of the atoms in a frame.
    #[arg(long = "unbuffered", default_value_t=true, action=clap::ArgAction::SetFalse)]
    is_buffered: bool,

    /// Write the trajectory in reverse.
    ///
    /// The direction of the selection is unaffected by this flag. To also reverse the frame
    /// selection, use `--reverse-frame-selection`.
    #[arg(long, short = 'r')]
    reverse: bool,

    /// Reverse the frame selection.
    ///
    /// Can be used independent of whether `--reverse` is set.
    #[arg(long, short = 'R')]
    reverse_frame_selection: bool,

    /// Print the time (ps) value for the selected frames to standard output.
    #[arg(long)]
    times: bool,

    /// Print the step number for the selected frames to standard output.
    ///
    /// If both `times` and `steps` are active, they will be separated by tabs and printed in that order.
    #[arg(long)]
    steps: bool,

    /// Force set the magic number of the output file.
    #[arg(long, hide = true)]
    force_magic: Option<i32>,
}

#[cfg(not(feature = "cli"))]
fn main() {
    panic!("The 'cli' feature must be enabled to use this binary.");
}

#[cfg(feature = "cli")]
fn main() -> std::io::Result<()> {
    let args = Args::parse();

    let file = std::fs::File::open(&args.input).unwrap_or_else(|err| {
        eprintln!(
            "ERROR: Failed to read trajectory from {:?}: {err}",
            &args.input
        );
        std::process::exit(1)
    });
    let mut reader = XTCReader::new(file);

    if args.info {
        let offsets = reader.determine_offsets(None)?;
        let name = args
            .input
            .file_name()
            .and_then(|name| name.to_str())
            .unwrap_or_default();
        println!("name:    {name}",);
        println!("path:    {:?}", &args.input);
        println!("nframes: {}", offsets.len());
        let headers = offsets
            .iter()
            .map(|&offset| -> io::Result<Header> {
                reader.file.seek(SeekFrom::Start(offset))?;
                reader.read_header()
            })
            .collect::<io::Result<Vec<_>>>()?;
        let natoms = headers
            .first()
            .map(|header| header.natoms.to_string())
            .unwrap_or("?".to_string());
        println!("natoms:  {natoms}");
        let (first, last) = (headers.first(), headers.last());

        let first_step = first.map(|header| header.step);
        let last_step = last.map(|header| header.step);
        let steps = match (first_step, last_step) {
            (None, None) => "?".to_string(),
            (None, Some(_)) => unreachable!(),
            (Some(first), None) => first.to_string(),
            (Some(first), Some(last)) => format!("{first}-{last}"),
        };
        println!("steps:   {steps}");

        let first_time = first.map(|header| header.time);
        let last_time = last.map(|header| header.time);
        let times = match (first_time, last_time) {
            (None, None) => "?".to_string(),
            (None, Some(_)) => unreachable!(),
            (Some(first), None) => first.to_string(),
            (Some(first), Some(last)) => format!("{first}-{last}"),
        };
        println!("time:    {times} ps");

        let magic = match first {
            Some(Header { magic, .. }) => magic.to_string(),
            None => "?".to_string(),
        };
        println!("magic:   {magic}");

        return Ok(());
    }

    let write = args
        .write
        .expect("write arguments must be available if --info is not passed");
    let mut writer = BufWriter::new(std::fs::File::create(&write.output).unwrap_or_else(|err| {
        eprintln!(
            "ERROR: Failed to write processed trajectory to {:?}: {err}",
            &write.output
        );
        std::process::exit(1)
    }));
    filter_frames(&mut reader, &mut writer, write)
}
