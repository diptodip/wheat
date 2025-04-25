#[derive(Copy, Clone)]
pub struct PRNG {
    pub state: u32,
}

impl PRNG {
    pub fn xor_shift32(&mut self) -> u32 {
        self.state ^= self.state << 13;
        self.state ^= self.state >> 17;
        self.state ^= self.state << 5;
        return self.state;
    }

    pub fn warm_up_xor_shift(&mut self) {
        for _ in 0..11 {
            let _x = self.xor_shift32();
        }
    }

    pub fn gen(&mut self) -> f32 {
        let x = self.xor_shift32();
        return x as f32 / u32::MAX as f32;
    }

    pub fn gen_range(&mut self, min: f32, max: f32) -> f32 {
        return self.gen() * (max - min) + min;
    }
}
