mod midi_gen;

use std::fs;
use std::path::PathBuf;

use clap::Parser;
use harp_voicing_wasm::abc;

const MODE_NAMES: [&str; 7] = [
    "ionian",
    "dorian",
    "phrygian",
    "lydian",
    "mixolydian",
    "aeolian",
    "locrian",
];

#[derive(Parser)]
#[command(name = "harp-midi", about = "Generate MIDI files from OpenHymnal ABC with harp voicings")]
struct Cli {
    /// Input ABC file path
    input: PathBuf,

    /// Output directory (default: ./output)
    #[arg(default_value = "output")]
    output_dir: PathBuf,

    /// Only generate hymn number N
    #[arg(long)]
    hymn: Option<String>,

    /// Only generate one mode (0=Ionian .. 6=Locrian)
    #[arg(long)]
    mode: Option<usize>,

    /// Override tempo (BPM)
    #[arg(long)]
    tempo: Option<u32>,

    /// MIDI velocity (0-127)
    #[arg(long, default_value = "80")]
    velocity: u8,
}

fn slugify(s: &str) -> String {
    s.to_lowercase()
        .chars()
        .map(|c| if c.is_alphanumeric() { c } else { '_' })
        .collect::<String>()
        .trim_matches('_')
        .to_string()
}

fn main() {
    let cli = Cli::parse();

    // Read as bytes first, then try UTF-8, falling back to Latin-1
    let raw_bytes = fs::read(&cli.input).unwrap_or_else(|e| {
        eprintln!("Error reading {}: {}", cli.input.display(), e);
        std::process::exit(1);
    });
    let abc_text = match String::from_utf8(raw_bytes.clone()) {
        Ok(s) => s,
        Err(_) => {
            // Latin-1: each byte maps directly to a Unicode code point
            raw_bytes.iter().map(|&b| b as char).collect()
        }
    };

    eprintln!("Parsing ABC file...");
    let mut hymns = abc::parse_lead_sheets(&abc_text);
    eprintln!("Found {} hymns", hymns.len());

    // Apply tempo override
    if let Some(bpm) = cli.tempo {
        for h in &mut hymns {
            h.tempo_bpm = bpm;
        }
    }

    // Filter by hymn number if specified
    if let Some(ref num) = cli.hymn {
        hymns.retain(|h| h.number == *num);
        if hymns.is_empty() {
            eprintln!("No hymn with number '{}' found", num);
            std::process::exit(1);
        }
    }

    // Determine which modes to generate
    let modes: Vec<usize> = match cli.mode {
        Some(m) if m < 7 => vec![m],
        Some(m) => {
            eprintln!("Invalid mode {}; must be 0-6", m);
            std::process::exit(1);
        }
        None => (0..7).collect(),
    };

    fs::create_dir_all(&cli.output_dir).unwrap_or_else(|e| {
        eprintln!("Error creating output dir: {}", e);
        std::process::exit(1);
    });

    let mut total = 0;
    for hymn in &hymns {
        let dir_name = format!(
            "{:0>3}_{}",
            hymn.number,
            slugify(&hymn.title).chars().take(30).collect::<String>()
        );
        let hymn_dir = cli.output_dir.join(&dir_name);
        fs::create_dir_all(&hymn_dir).ok();

        for &mode in &modes {
            let midi_data = midi_gen::hymn_to_midi(hymn, mode as i16, cli.velocity);
            let filename = format!(
                "{:0>3}_{}.mid",
                hymn.number,
                MODE_NAMES[mode]
            );
            let path = hymn_dir.join(&filename);
            fs::write(&path, &midi_data).unwrap_or_else(|e| {
                eprintln!("Error writing {}: {}", path.display(), e);
            });
            total += 1;
        }

        eprint!("\r{}: {} ({} files)", hymn.number, hymn.title, total);
    }

    eprintln!("\nDone. Generated {} MIDI files.", total);
}
