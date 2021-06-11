pub fn output_ppm(image: Vec<Vec<f64>>, rows: usize, cols: usize) {
    println!("P3");
    print!("{}", cols);
    print!(" ");
    println!("{}", rows);
    println!("{}", 255);
    let cap = 255.999;
    for i in 0..rows*cols {
        let row = i / cols;
        let col = i % rows;
        let r = (cap * image[i][0].sqrt().min(0.999).max(0.0)) as i32;
        let g = (cap * image[i][1].sqrt().min(0.999).max(0.0)) as i32;
        let b = (cap * image[i][2].sqrt().min(0.999).max(0.0)) as i32;
        println!("{} {} {} ", r, g, b);
    }
}
