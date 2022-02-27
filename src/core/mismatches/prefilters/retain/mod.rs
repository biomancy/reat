use std::ops::Range;

use bio_types::genome::{Interval, Position};

pub use rois::RetainROIFromList;
pub use sites::RetainSitesFromIntervals;

mod rois;
mod sites;

pub trait SitesRetainer {
    fn filter(&self, contig: &str, range: Range<Position>) -> Vec<Range<Position>>;
}

pub trait ROIRetainer {
    fn filter(&self, interval: &Interval, name: &str) -> bool;
}

// pub trait RetainableROIs {
//     fn filter(&self, contig: &str, range: Range<Position>) -> bool;
// }

// Имеет ли смысл тогда вообще делать хоть какие-то интерфейсы? Ведь один из них будет возвращать одно, а другой другое.
// Смысла в интерфейсе нет никакого.

// Это фильтры по локации, верно? Я не вижу никаких других вариантов.
// Ок, тогда там будет простая чушь, которая по одному чекает сайты/интервалы и потом их как-то обрабатывает.
//
//
//
// Как разбивать маленькие фигни? Через хеш? Или interval tree? Или еще как-то?
// Наверное, хеш будет быстрее. Или нет? Черт его знает.
// Вероятно да, быстрее.
