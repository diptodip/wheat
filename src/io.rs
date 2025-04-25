use std::fs::File;
use std::io::{BufWriter, Error, Write};

pub fn write_ppm(image: Vec<Vec<f32>>, rows: usize, cols: usize) -> Result<(), Error> {
    let mut f = BufWriter::new(File::create("image.ppm")?);
    writeln!(f, "P3")?;
    write!(f, "{}", cols)?;
    write!(f, " ")?;
    writeln!(f, "{}", rows)?;
    writeln!(f, "{}", 255)?;
    let cap = 255.999;
    for i in 0..rows * cols {
        let r = (cap * image[i][0].sqrt().min(0.999).max(0.0)) as i32;
        let g = (cap * image[i][1].sqrt().min(0.999).max(0.0)) as i32;
        let b = (cap * image[i][2].sqrt().min(0.999).max(0.0)) as i32;
        writeln!(f, "{} {} {} ", r, g, b)?;
    }
    f.flush()?;
    Ok(())
}
