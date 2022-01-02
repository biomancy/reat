pub trait Table {
    const LENGTH: usize;

    fn row(&self) -> Vec<String>;
    fn header() -> Vec<String>;
}
