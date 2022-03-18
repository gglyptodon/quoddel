fn main() {
    if let Err(e) = quoddel::get_args().and_then(quoddel::run) {
        eprintln!("{}", e);
    }
}
