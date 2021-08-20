use std::cell::RefCell;
use std::sync::Mutex;

use thread_local::ThreadLocal;

pub struct ThreadCache<Builder, Type>
where
    Builder: Fn() -> RefCell<Type>,
    Type: Send,
{
    mutex: Mutex<Builder>,
    thrlocal: ThreadLocal<RefCell<Type>>,
}

impl<Builder, Type> ThreadCache<Builder, Type>
where
    Builder: Fn() -> RefCell<Type>,
    Type: Send,
{
    pub fn new(builder: Builder) -> Self {
        Self { mutex: Mutex::new(builder), thrlocal: Default::default() }
    }

    pub fn get(&self) -> &RefCell<Type> {
        self.thrlocal.get_or(|| self.mutex.lock().unwrap()())
    }

    pub fn content(self) -> <ThreadLocal<RefCell<Type>> as IntoIterator>::IntoIter {
        self.thrlocal.into_iter()
    }
}
