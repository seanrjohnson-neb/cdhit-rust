fn main() {
    let data = std::fs::read(std::env::args().nth(1).unwrap()).unwrap();
    let n: i32 = std::env::args().nth(2).unwrap().parse().unwrap();
    let out = std::env::args().nth(3).unwrap();
    for (i, seg) in cdhit_core::divide(&data, n).into_iter().enumerate() {
        std::fs::write(format!("{out}-{i}"), seg).unwrap();
    }
}
