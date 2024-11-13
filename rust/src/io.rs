extern crate atty;
use std::collections::HashSet;
use std::ffi::OsStr;
use std::fs::{create_dir_all, File};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Result, Seek, Write};
use std::path::{Path, PathBuf};

use csv::ReaderBuilder;
use flate2::read::GzDecoder;
use flate2::write;
use flate2::Compression;

fn read_stdin() -> Vec<Vec<u8>> {
    let stdin = io::stdin();
    let mut list: Vec<Vec<u8>> = vec![];
    if atty::is(atty::Stream::Stdin) {
        eprintln!("No input on STDIN!");
        return list;
    }
    for line in stdin.lock().lines() {
        let line_as_vec = match line {
            Err(why) => panic!("couldn't read line: {}", why),
            Ok(l) => l.as_bytes().to_vec(),
        };
        list.push(line_as_vec)
    }
    list
}

pub fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where
    P: AsRef<Path>,
{
    let file = File::open(filename).expect("no such file");
    Ok(io::BufReader::new(file).lines())
}

fn read_file(file_path: &PathBuf) -> Vec<Vec<u8>> {
    let mut output: Vec<Vec<u8>> = vec![];
    if let Ok(lines) = read_lines(file_path) {
        for line in lines {
            let line_as_vec = match line {
                Err(why) => panic!("couldn't read line: {}", why),
                Ok(l) => l.as_bytes().to_vec(),
            };
            output.push(line_as_vec)
        }
    }
    output
}

pub fn get_list(file_path: &Option<PathBuf>) -> HashSet<Vec<u8>> {
    let list = match file_path {
        None => vec![],
        Some(p) if p == Path::new("-") => read_stdin(),
        Some(_) => read_file(file_path.as_ref().unwrap()),
    };
    HashSet::from_iter(list)
}

pub fn get_file_writer(file_path: &PathBuf) -> Box<dyn Write> {
    let file = match File::create(file_path) {
        Err(why) => panic!("couldn't open {}: {}", file_path.display(), why),
        Ok(file) => file,
    };

    let writer: Box<dyn Write> = if file_path.extension() == Some(OsStr::new("gz")) {
        Box::new(BufWriter::with_capacity(
            128 * 1024,
            write::GzEncoder::new(file, Compression::default()),
        ))
    } else {
        Box::new(BufWriter::with_capacity(128 * 1024, file))
    };
    writer
}

pub fn get_writer(file_path: &Option<PathBuf>) -> Box<dyn Write> {
    let writer: Box<dyn Write> = match file_path {
        Some(path) if path == Path::new("-") => Box::new(BufWriter::new(io::stdout().lock())),
        Some(path) => {
            create_dir_all(path.parent().unwrap()).unwrap();
            get_file_writer(path)
        }
        None => Box::new(BufWriter::new(io::stdout().lock())),
    };
    writer
}

pub fn write_list(entries: &HashSet<Vec<u8>>, file_path: &Option<PathBuf>) -> Result<()> {
    let mut writer = get_writer(file_path);
    for line in entries.iter() {
        writeln!(&mut writer, "{}", String::from_utf8(line.to_vec()).unwrap())?;
    }
    Ok(())
}

pub fn append_to_path(p: &PathBuf, s: &str) -> PathBuf {
    let mut p = p.clone().into_os_string();
    p.push(s);
    p.into()
}

pub fn file_reader(path: PathBuf) -> Option<Box<dyn BufRead>> {
    let file = File::open(&path).expect("no such file");

    if path.ends_with(".gz") {
        return Some(Box::new(BufReader::new(GzDecoder::new(file))));
    } else {
        return Some(Box::new(BufReader::new(file)));
    };
}

pub fn csv_reader(
    header: bool,
    delimiter: u8,
    path: PathBuf,
) -> Result<csv::Reader<Box<dyn Read>>> {
    let file = File::open(&path)?;
    let mut buf = [0; 2];
    let mut reader = file.try_clone()?;

    reader.read_exact(&mut buf)?;
    reader.rewind()?;

    let rdr: csv::Reader<Box<dyn std::io::Read>> = if buf == [0x1f, 0x8b] {
        ReaderBuilder::new()
            .has_headers(header)
            .delimiter(delimiter)
            .from_reader(Box::new(BufReader::new(GzDecoder::new(file))))
    } else {
        ReaderBuilder::new()
            .has_headers(header)
            .delimiter(delimiter)
            .from_reader(Box::new(BufReader::new(file)))
    };

    Ok(rdr)
}
