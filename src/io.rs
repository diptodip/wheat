pub fn output_ppm(image: Vec<Vec<Vec<f64>>>, rows: usize, cols: usize) {
    let rows = rows;
    let cols = cols;
    println!("P3");
    print!("{}", cols);
    print!(" ");
    println!("{}", rows);
    println!("{}", 255);
    let cap = 255.999;
    for i in 0..rows {
        for j in 0..cols {
            let r = (cap * image[i][j][0].sqrt().min(0.999).max(0.0)) as i32;
            let g = (cap * image[i][j][1].sqrt().min(0.999).max(0.0)) as i32;
            let b = (cap * image[i][j][2].sqrt().min(0.999).max(0.0)) as i32;
            println!("{} {} {} ", r, g, b);
        }
    }
}
