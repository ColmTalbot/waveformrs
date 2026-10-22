// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL and fdlibm: see NOTICE.md for their authors.
//! A minimal, read-only reader for the subset of HDF5 used by the
//! LALSimulation surrogate data files (e.g. `NRSur7dq4_v1.0.h5`).
//!
//! The files are written with h5py's `libver="earliest"` layout, which uses
//! the original (HDF5 1.0) on-disk structures only. This module supports
//! exactly that subset:
//!
//! * superblock versions 0 and 1 (also behind a user block);
//! * version 1 object headers, including continuation blocks;
//! * "old-style" groups: a symbol table message pointing at a version 1
//!   B-tree of symbol table nodes (`SNOD`) and a local heap for the names;
//! * datasets with a version 1 or 2 dataspace (scalar or simple), an IEEE
//!   floating-point (4 or 8 byte) or integer (1, 2, 4 or 8 byte) datatype in
//!   either byte order, and a version 3 contiguous data layout without
//!   filters.
//!
//! Anything else (version 2+ superblocks or object headers, compact or
//! indexed "new-style" groups, chunked, compact, external or filtered
//! (compressed) storage, shared or committed datatypes, compound, string or
//! variable-length types) is rejected with [`Hdf5Error::Unsupported`]
//! rather than misread. Attributes are ignored.
//!
//! The whole file is read into memory once ([`File::open`]); groups and
//! datasets are views into that buffer, so opening many datasets is cheap.

use std::fmt;
use std::path::Path;

/// Errors from reading an HDF5 file.
#[derive(Debug)]
pub enum Hdf5Error {
    /// The file could not be read.
    Io(std::io::Error),
    /// The file is not valid HDF5 (or is truncated/corrupt).
    Format(String),
    /// A valid HDF5 feature outside the supported subset.
    Unsupported(String),
    /// No group member of that name.
    NotFound(String),
    /// The object exists but has a different kind, element type or shape
    /// than requested.
    Type(String),
}

impl fmt::Display for Hdf5Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(err) => write!(f, "HDF5 I/O error: {err}"),
            Self::Format(why) => write!(f, "invalid HDF5 file: {why}"),
            Self::Unsupported(what) => write!(f, "unsupported HDF5 feature: {what}"),
            Self::NotFound(name) => write!(f, "HDF5 object not found: {name}"),
            Self::Type(why) => write!(f, "HDF5 type mismatch: {why}"),
        }
    }
}

impl std::error::Error for Hdf5Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io(err) => Some(err),
            _ => None,
        }
    }
}

impl From<std::io::Error> for Hdf5Error {
    fn from(err: std::io::Error) -> Self {
        Self::Io(err)
    }
}

pub type Result<T> = std::result::Result<T, Hdf5Error>;

fn format_error<T>(why: impl Into<String>) -> Result<T> {
    Err(Hdf5Error::Format(why.into()))
}

fn unsupported<T>(what: impl Into<String>) -> Result<T> {
    Err(Hdf5Error::Unsupported(what.into()))
}

const SIGNATURE: [u8; 8] = [0x89, b'H', b'D', b'F', b'\r', b'\n', 0x1a, b'\n'];

// Object header message types.
const MSG_NIL: u16 = 0x0000;
const MSG_DATASPACE: u16 = 0x0001;
const MSG_LINK_INFO: u16 = 0x0002;
const MSG_DATATYPE: u16 = 0x0003;
const MSG_EXTERNAL_FILES: u16 = 0x0007;
const MSG_LAYOUT: u16 = 0x0008;
const MSG_GROUP_INFO: u16 = 0x000A;
const MSG_FILTER_PIPELINE: u16 = 0x000B;
const MSG_LINK: u16 = 0x0006;
const MSG_CONTINUATION: u16 = 0x0010;
const MSG_SYMBOL_TABLE: u16 = 0x0011;

/// An HDF5 file held in memory.
pub struct File {
    bytes: Vec<u8>,
    /// Absolute position of the base address (start of the superblock).
    base: usize,
    size_of_offsets: usize,
    size_of_lengths: usize,
    root_header: u64,
}

impl fmt::Debug for File {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("File")
            .field("len", &self.bytes.len())
            .field("base", &self.base)
            .finish()
    }
}

/// Bounds-checked little-endian reads from the file buffer.
#[derive(Clone, Copy)]
struct Cursor<'a> {
    bytes: &'a [u8],
    pos: usize,
}

impl<'a> Cursor<'a> {
    fn take(&mut self, n: usize) -> Result<&'a [u8]> {
        let end = self
            .pos
            .checked_add(n)
            .filter(|&end| end <= self.bytes.len());
        match end {
            Some(end) => {
                let slice = &self.bytes[self.pos..end];
                self.pos = end;
                Ok(slice)
            }
            None => format_error(format!(
                "read of {n} bytes at {} beyond the end of the data ({})",
                self.pos,
                self.bytes.len()
            )),
        }
    }

    fn u8(&mut self) -> Result<u8> {
        Ok(self.take(1)?[0])
    }

    fn u16(&mut self) -> Result<u16> {
        Ok(u16::from_le_bytes(self.take(2)?.try_into().unwrap()))
    }

    fn u32(&mut self) -> Result<u32> {
        Ok(u32::from_le_bytes(self.take(4)?.try_into().unwrap()))
    }

    /// An unsigned integer of `n` bytes (1 to 8).
    fn uint(&mut self, n: usize) -> Result<u64> {
        let bytes = self.take(n)?;
        let mut value = 0u64;
        for (i, &b) in bytes.iter().enumerate() {
            value |= (b as u64) << (8 * i);
        }
        Ok(value)
    }

    fn skip(&mut self, n: usize) -> Result<()> {
        self.take(n).map(|_| ())
    }
}

/// The "undefined address" (all bits set) for `size` bytes.
fn undefined_address(size: usize) -> u64 {
    if size >= 8 {
        u64::MAX
    } else {
        (1u64 << (8 * size)) - 1
    }
}

impl File {
    /// Read the whole file and parse its superblock.
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        Self::from_bytes(std::fs::read(path)?)
    }

    /// Parse an HDF5 file held in memory.
    pub fn from_bytes(bytes: Vec<u8>) -> Result<Self> {
        // The superblock is at 0, 512, 1024, 2048, ... (after a user block).
        let mut location = 0usize;
        loop {
            if location + SIGNATURE.len() > bytes.len() {
                return format_error("HDF5 signature not found");
            }
            if bytes[location..location + SIGNATURE.len()] == SIGNATURE {
                break;
            }
            location = if location == 0 { 512 } else { 2 * location };
        }
        let mut cursor = Cursor {
            bytes: &bytes,
            pos: location + SIGNATURE.len(),
        };
        let version = cursor.u8()?;
        if version > 1 {
            return unsupported(format!(
                "superblock version {version} (only versions 0 and 1 are supported)"
            ));
        }
        let free_space_version = cursor.u8()?;
        let root_entry_version = cursor.u8()?;
        cursor.skip(1)?;
        let shared_header_version = cursor.u8()?;
        if free_space_version != 0 || root_entry_version != 0 || shared_header_version != 0 {
            return unsupported(format!(
                "superblock sub-versions {free_space_version}/{root_entry_version}/\
                 {shared_header_version}"
            ));
        }
        let size_of_offsets = cursor.u8()? as usize;
        let size_of_lengths = cursor.u8()? as usize;
        for size in [size_of_offsets, size_of_lengths] {
            if !matches!(size, 2 | 4 | 8) {
                return format_error(format!("invalid size of offsets/lengths {size}"));
            }
        }
        cursor.skip(1)?;
        let _group_leaf_k = cursor.u16()?;
        let _group_internal_k = cursor.u16()?;
        let _consistency_flags = cursor.u32()?;
        if version == 1 {
            let _indexed_storage_k = cursor.u16()?;
            cursor.skip(2)?;
        }
        let base_address = cursor.uint(size_of_offsets)?;
        let _free_space_address = cursor.uint(size_of_offsets)?;
        let _end_of_file = cursor.uint(size_of_offsets)?;
        let _driver_info = cursor.uint(size_of_offsets)?;
        // root group symbol table entry
        let _link_name_offset = cursor.uint(size_of_offsets)?;
        let root_header = cursor.uint(size_of_offsets)?;
        let base = usize::try_from(base_address)
            .ok()
            .filter(|&b| b <= bytes.len())
            .ok_or_else(|| Hdf5Error::Format(format!("invalid base address {base_address}")))?;
        Ok(Self {
            bytes,
            base,
            size_of_offsets,
            size_of_lengths,
            root_header,
        })
    }

    /// A cursor at file address `address` (relative to the base address).
    fn at(&self, address: u64) -> Result<Cursor<'_>> {
        let pos = usize::try_from(address)
            .ok()
            .and_then(|a| a.checked_add(self.base))
            .filter(|&p| p <= self.bytes.len());
        match pos {
            Some(pos) => Ok(Cursor {
                bytes: &self.bytes,
                pos,
            }),
            None => format_error(format!("address {address} outside the file")),
        }
    }

    fn offset(&self, cursor: &mut Cursor) -> Result<u64> {
        cursor.uint(self.size_of_offsets)
    }

    fn length(&self, cursor: &mut Cursor) -> Result<u64> {
        cursor.uint(self.size_of_lengths)
    }

    /// The root group.
    pub fn root(&self) -> Result<Group<'_>> {
        Group::open(self, "/".to_string(), self.root_header)
    }

    /// The group at a `/`-separated path from the root.
    pub fn group(&self, path: &str) -> Result<Group<'_>> {
        let mut group = self.root()?;
        for part in path.split('/').filter(|p| !p.is_empty()) {
            group = group.group(part)?;
        }
        Ok(group)
    }

    /// The dataset at a `/`-separated path from the root.
    pub fn dataset(&self, path: &str) -> Result<Dataset<'_>> {
        let (dir, name) = match path.trim_start_matches('/').rsplit_once('/') {
            Some((dir, name)) => (dir, name),
            None => ("", path.trim_start_matches('/')),
        };
        self.group(dir)?.dataset(name)
    }

    /// Parse the messages of the version 1 object header at `address`.
    fn object_header(&self, address: u64) -> Result<Vec<Message<'_>>> {
        let mut cursor = self.at(address)?;
        if cursor.bytes.get(cursor.pos..cursor.pos + 4) == Some(b"OHDR") {
            return unsupported("version 2 object headers");
        }
        let version = cursor.u8()?;
        if version != 1 {
            return format_error(format!(
                "object header version {version} at address {address}"
            ));
        }
        cursor.skip(1)?;
        let n_messages = cursor.u16()? as usize;
        let _reference_count = cursor.u32()?;
        let header_size = cursor.u32()? as u64;
        cursor.skip(4)?; // padding to an 8-byte boundary
        let mut blocks = vec![(address + 16, header_size)];
        let mut messages = Vec::with_capacity(n_messages);
        let mut block_index = 0;
        while block_index < blocks.len() {
            let (start, size) = blocks[block_index];
            block_index += 1;
            let mut cursor = self.at(start)?;
            let end = cursor
                .pos
                .checked_add(size as usize)
                .filter(|&e| e <= self.bytes.len())
                .ok_or_else(|| Hdf5Error::Format("object header block beyond EOF".into()))?;
            while cursor.pos + 8 <= end {
                let kind = cursor.u16()?;
                let data_size = cursor.u16()? as usize;
                let flags = cursor.u8()?;
                cursor.skip(3)?;
                if cursor.pos + data_size > end {
                    return format_error("object header message overruns its block");
                }
                let data = &self.bytes[cursor.pos..cursor.pos + data_size];
                cursor.skip(data_size)?;
                if kind == MSG_NIL {
                    continue;
                }
                if flags & 0x02 != 0 {
                    return unsupported(format!("shared object header message (type {kind:#x})"));
                }
                if kind == MSG_CONTINUATION {
                    let mut c = Cursor {
                        bytes: data,
                        pos: 0,
                    };
                    let offset = self.offset(&mut c)?;
                    let length = self.length(&mut c)?;
                    blocks.push((offset, length));
                    continue;
                }
                messages.push(Message { kind, data });
            }
        }
        Ok(messages)
    }
}

struct Message<'a> {
    kind: u16,
    data: &'a [u8],
}

/// Kind of a group member.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ObjectKind {
    Group,
    Dataset,
}

/// An (old-style) HDF5 group: its members in name order.
#[derive(Debug)]
pub struct Group<'f> {
    file: &'f File,
    path: String,
    /// `(name, object header address)` sorted by name.
    entries: Vec<(String, u64)>,
}

impl<'f> Group<'f> {
    fn open(file: &'f File, path: String, header: u64) -> Result<Self> {
        let messages = file.object_header(header)?;
        let mut symbol_table = None;
        for message in &messages {
            match message.kind {
                MSG_SYMBOL_TABLE => {
                    let mut c = Cursor {
                        bytes: message.data,
                        pos: 0,
                    };
                    let btree = file.offset(&mut c)?;
                    let heap = file.offset(&mut c)?;
                    symbol_table = Some((btree, heap));
                }
                MSG_LINK | MSG_LINK_INFO | MSG_GROUP_INFO => {
                    return unsupported(format!("new-style (link message) group {path}"));
                }
                _ => {}
            }
        }
        let Some((btree, heap)) = symbol_table else {
            return Err(Hdf5Error::Type(format!("{path} is not a group")));
        };
        let heap_data = Self::local_heap(file, heap)?;
        let mut entries = Vec::new();
        Self::collect(file, btree, heap_data, &mut entries, 0)?;
        entries.sort_by(|a, b| a.0.cmp(&b.0));
        Ok(Self {
            file,
            path,
            entries,
        })
    }

    /// The data segment of a local heap.
    fn local_heap(file: &File, address: u64) -> Result<&[u8]> {
        let mut c = file.at(address)?;
        if c.take(4)? != b"HEAP" {
            return format_error(format!("no local heap at address {address}"));
        }
        let version = c.u8()?;
        if version != 0 {
            return unsupported(format!("local heap version {version}"));
        }
        c.skip(3)?;
        let size = file.length(&mut c)? as usize;
        let _free_list = file.length(&mut c)?;
        let data_address = file.offset(&mut c)?;
        let mut data = file.at(data_address)?;
        data.take(size)
    }

    /// Walk a version 1 group B-tree, collecting the symbol table entries.
    fn collect(
        file: &File,
        address: u64,
        heap: &[u8],
        entries: &mut Vec<(String, u64)>,
        depth: usize,
    ) -> Result<()> {
        if depth > 64 {
            return format_error("B-tree too deep (cycle?)");
        }
        let mut c = file.at(address)?;
        if c.take(4)? != b"TREE" {
            return format_error(format!("no B-tree node at address {address}"));
        }
        let node_type = c.u8()?;
        if node_type != 0 {
            return format_error(format!("B-tree node type {node_type} in a group"));
        }
        let level = c.u8()?;
        let n_entries = c.u16()? as usize;
        let _left = file.offset(&mut c)?;
        let _right = file.offset(&mut c)?;
        for _ in 0..n_entries {
            let _key = file.length(&mut c)?;
            let child = file.offset(&mut c)?;
            if level > 0 {
                Self::collect(file, child, heap, entries, depth + 1)?;
            } else {
                Self::symbol_node(file, child, heap, entries)?;
            }
        }
        Ok(())
    }

    fn symbol_node(
        file: &File,
        address: u64,
        heap: &[u8],
        entries: &mut Vec<(String, u64)>,
    ) -> Result<()> {
        let mut c = file.at(address)?;
        if c.take(4)? != b"SNOD" {
            return format_error(format!("no symbol table node at address {address}"));
        }
        let version = c.u8()?;
        if version != 1 {
            return unsupported(format!("symbol table node version {version}"));
        }
        c.skip(1)?;
        let n_symbols = c.u16()? as usize;
        for _ in 0..n_symbols {
            let name_offset = file.length(&mut c)? as usize;
            let header = file.offset(&mut c)?;
            let _cache_type = c.u32()?;
            c.skip(4 + 16)?;
            let name = heap
                .get(name_offset..)
                .and_then(|s| s.iter().position(|&b| b == 0).map(|end| &s[..end]))
                .ok_or_else(|| Hdf5Error::Format("link name outside the local heap".into()))?;
            let name = String::from_utf8(name.to_vec())
                .map_err(|_| Hdf5Error::Format("link name is not UTF-8".into()))?;
            entries.push((name, header));
        }
        Ok(())
    }

    /// The path of this group.
    pub fn path(&self) -> &str {
        &self.path
    }

    /// The names of the members, sorted.
    pub fn member_names(&self) -> impl Iterator<Item = &str> {
        self.entries.iter().map(|(name, _)| name.as_str())
    }

    /// The number of members.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Whether the group has no members.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    fn child_path(&self, name: &str) -> String {
        if self.path == "/" {
            format!("/{name}")
        } else {
            format!("{}/{name}", self.path)
        }
    }

    fn header_of(&self, name: &str) -> Result<u64> {
        self.entries
            .binary_search_by(|(n, _)| n.as_str().cmp(name))
            .map(|i| self.entries[i].1)
            .map_err(|_| Hdf5Error::NotFound(self.child_path(name)))
    }

    /// Whether a member is a group or a dataset.
    pub fn kind(&self, name: &str) -> Result<ObjectKind> {
        let messages = self.file.object_header(self.header_of(name)?)?;
        if messages.iter().any(|m| m.kind == MSG_SYMBOL_TABLE) {
            Ok(ObjectKind::Group)
        } else if messages.iter().any(|m| m.kind == MSG_LAYOUT) {
            Ok(ObjectKind::Dataset)
        } else if messages
            .iter()
            .any(|m| matches!(m.kind, MSG_LINK | MSG_LINK_INFO | MSG_GROUP_INFO))
        {
            unsupported(format!("new-style group {}", self.child_path(name)))
        } else {
            unsupported(format!(
                "object {} (committed datatype?)",
                self.child_path(name)
            ))
        }
    }

    /// The sub-group `name`.
    pub fn group(&self, name: &str) -> Result<Group<'f>> {
        Group::open(self.file, self.child_path(name), self.header_of(name)?)
    }

    /// The dataset `name`.
    pub fn dataset(&self, name: &str) -> Result<Dataset<'f>> {
        Dataset::open(self.file, self.child_path(name), self.header_of(name)?)
    }
}

/// Element type of a dataset.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Datatype {
    /// IEEE 754 binary floating point of `size` (4 or 8) bytes.
    Float { size: usize, big_endian: bool },
    /// Two's complement or unsigned integer of `size` (1, 2, 4 or 8) bytes.
    Integer {
        size: usize,
        signed: bool,
        big_endian: bool,
    },
}

impl Datatype {
    fn size(&self) -> usize {
        match *self {
            Self::Float { size, .. } | Self::Integer { size, .. } => size,
        }
    }

    fn parse(data: &[u8]) -> Result<Self> {
        let mut c = Cursor {
            bytes: data,
            pos: 0,
        };
        let class_version = c.u8()?;
        let class = class_version & 0x0f;
        let version = class_version >> 4;
        let bits = c.uint(3)? as u32;
        let size = c.u32()? as usize;
        if version != 1 {
            return unsupported(format!("datatype message version {version}"));
        }
        match class {
            0 => {
                let big_endian = bits & 0x01 != 0;
                let signed = bits & 0x08 != 0;
                let offset = c.u16()?;
                let precision = c.u16()? as usize;
                if !matches!(size, 1 | 2 | 4 | 8) || offset != 0 || precision != 8 * size {
                    return unsupported(format!(
                        "integer type of {size} bytes with bit offset {offset} and precision \
                         {precision}"
                    ));
                }
                Ok(Self::Integer {
                    size,
                    signed,
                    big_endian,
                })
            }
            1 => {
                if bits & 0x40 != 0 {
                    return unsupported("VAX floating-point byte order");
                }
                let big_endian = bits & 0x01 != 0;
                let normalization = (bits >> 4) & 0x03;
                let sign = (bits >> 8) & 0xff;
                let offset = c.u16()?;
                let precision = c.u16()? as u32;
                let exponent_location = c.u8()? as u32;
                let exponent_size = c.u8()? as u32;
                let mantissa_location = c.u8()? as u32;
                let mantissa_size = c.u8()? as u32;
                let bias = c.u32()?;
                let ieee = match size {
                    8 => (64, 52, 11, 0, 52, 1023, 63),
                    4 => (32, 23, 8, 0, 23, 127, 31),
                    _ => (0, 0, 0, 0, 0, 0, 0),
                };
                let layout = (
                    precision,
                    exponent_location,
                    exponent_size,
                    mantissa_location,
                    mantissa_size,
                    bias,
                    sign,
                );
                if offset != 0 || layout != ieee || normalization != 2 {
                    return unsupported(format!(
                        "non-IEEE floating-point type of {size} bytes {layout:?}"
                    ));
                }
                Ok(Self::Float { size, big_endian })
            }
            _ => unsupported(format!("datatype class {class}")),
        }
    }
}

/// A dataset with contiguous storage.
#[derive(Debug)]
pub struct Dataset<'f> {
    file: &'f File,
    path: String,
    shape: Vec<usize>,
    dtype: Datatype,
    /// Absolute byte range of the data in the file buffer.
    data: std::ops::Range<usize>,
}

impl<'f> Dataset<'f> {
    fn open(file: &'f File, path: String, header: u64) -> Result<Self> {
        let messages = file.object_header(header)?;
        let mut shape = None;
        let mut dtype = None;
        let mut layout = None;
        for message in &messages {
            let mut c = Cursor {
                bytes: message.data,
                pos: 0,
            };
            match message.kind {
                MSG_DATASPACE => {
                    let version = c.u8()?;
                    let rank = c.u8()? as usize;
                    let flags = c.u8()?;
                    let kind = match version {
                        1 => {
                            c.skip(5)?;
                            if rank == 0 {
                                0
                            } else {
                                1
                            }
                        }
                        2 => c.u8()?,
                        _ => return unsupported(format!("dataspace version {version}")),
                    };
                    if kind == 2 {
                        return unsupported(format!("null dataspace in {path}"));
                    }
                    let _ = flags;
                    let mut dims = Vec::with_capacity(rank);
                    for _ in 0..rank {
                        let dim = file.length(&mut c)?;
                        dims.push(usize::try_from(dim).map_err(|_| {
                            Hdf5Error::Format(format!("dimension {dim} too large"))
                        })?);
                    }
                    shape = Some(dims);
                }
                MSG_DATATYPE => dtype = Some(Datatype::parse(message.data)?),
                MSG_LAYOUT => {
                    let version = c.u8()?;
                    if version != 3 {
                        return unsupported(format!(
                            "data layout message version {version} in {path}"
                        ));
                    }
                    let class = c.u8()?;
                    match class {
                        1 => {
                            let address = file.offset(&mut c)?;
                            let size = file.length(&mut c)?;
                            layout = Some((address, size));
                        }
                        0 => return unsupported(format!("compact storage in {path}")),
                        2 => return unsupported(format!("chunked storage in {path}")),
                        _ => return unsupported(format!("layout class {class} in {path}")),
                    }
                }
                MSG_FILTER_PIPELINE => {
                    return unsupported(format!("filtered (compressed) data in {path}"))
                }
                MSG_EXTERNAL_FILES => return unsupported(format!("external storage in {path}")),
                MSG_SYMBOL_TABLE => return Err(Hdf5Error::Type(format!("{path} is a group"))),
                _ => {}
            }
        }
        let (Some(shape), Some(dtype), Some((address, size))) = (shape, dtype, layout) else {
            return Err(Hdf5Error::Type(format!(
                "{path} is not a dataset (missing dataspace, datatype or layout)"
            )));
        };
        let n_elements = shape
            .iter()
            .try_fold(1usize, |n, &d| n.checked_mul(d))
            .ok_or_else(|| Hdf5Error::Format(format!("{path}: too many elements")))?;
        let n_bytes = n_elements * dtype.size();
        let data = if n_bytes == 0 {
            0..0
        } else {
            if address == undefined_address(file.size_of_offsets) {
                return format_error(format!("{path}: data was never written"));
            }
            if size < n_bytes as u64 {
                return format_error(format!(
                    "{path}: storage of {size} bytes for {n_bytes} bytes of data"
                ));
            }
            let start = file.at(address)?.pos;
            if start + n_bytes > file.bytes.len() {
                return format_error(format!("{path}: data beyond the end of the file"));
            }
            start..start + n_bytes
        };
        Ok(Self {
            file,
            path,
            shape,
            dtype,
            data,
        })
    }

    /// The full path of the dataset.
    pub fn path(&self) -> &str {
        &self.path
    }

    /// The dimensions (empty for a scalar).
    pub fn shape(&self) -> &[usize] {
        &self.shape
    }

    /// The number of elements.
    pub fn len(&self) -> usize {
        self.shape.iter().product()
    }

    /// Whether the dataset has no elements.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The element type.
    pub fn dtype(&self) -> Datatype {
        self.dtype
    }

    /// The raw bytes of the data (in the file's byte order).
    pub fn raw(&self) -> &'f [u8] {
        &self.file.bytes[self.data.clone()]
    }

    /// The data as `f64` (row-major), from a floating-point dataset.
    pub fn read_f64(&self) -> Result<Vec<f64>> {
        let raw = self.raw();
        match self.dtype {
            Datatype::Float {
                size: 8,
                big_endian,
            } => Ok(raw
                .chunks_exact(8)
                .map(|b| {
                    let b: [u8; 8] = b.try_into().unwrap();
                    if big_endian {
                        f64::from_be_bytes(b)
                    } else {
                        f64::from_le_bytes(b)
                    }
                })
                .collect()),
            Datatype::Float {
                size: 4,
                big_endian,
            } => Ok(raw
                .chunks_exact(4)
                .map(|b| {
                    let b: [u8; 4] = b.try_into().unwrap();
                    (if big_endian {
                        f32::from_be_bytes(b)
                    } else {
                        f32::from_le_bytes(b)
                    }) as f64
                })
                .collect()),
            other => Err(Hdf5Error::Type(format!(
                "{}: requested f64 from {other:?}",
                self.path
            ))),
        }
    }

    /// The data as `i64` (row-major), from an integer dataset.
    pub fn read_i64(&self) -> Result<Vec<i64>> {
        let Datatype::Integer {
            size,
            signed,
            big_endian,
        } = self.dtype
        else {
            return Err(Hdf5Error::Type(format!(
                "{}: requested i64 from {:?}",
                self.path, self.dtype
            )));
        };
        let raw = self.raw();
        if size == 8 && signed && !big_endian {
            // the common case, without per-element branches
            return Ok(raw
                .chunks_exact(8)
                .map(|b| i64::from_le_bytes(b.try_into().unwrap()))
                .collect());
        }
        let mut out = Vec::with_capacity(raw.len() / size);
        for b in raw.chunks_exact(size) {
            let mut le = [0u8; 8];
            if big_endian {
                for (dst, src) in le.iter_mut().zip(b.iter().rev()) {
                    *dst = *src;
                }
            } else {
                le[..size].copy_from_slice(b);
            }
            if signed && le[size - 1] & 0x80 != 0 {
                le[size..].fill(0xff);
            }
            let value = u64::from_le_bytes(le);
            if !signed && value > i64::MAX as u64 {
                return Err(Hdf5Error::Type(format!(
                    "{}: unsigned value {value} does not fit in i64",
                    self.path
                )));
            }
            out.push(value as i64);
        }
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_non_hdf5() {
        assert!(matches!(
            File::from_bytes(b"not an hdf5 file at all".to_vec()),
            Err(Hdf5Error::Format(_))
        ));
    }

    #[test]
    fn rejects_new_superblocks() {
        let mut bytes = SIGNATURE.to_vec();
        bytes.extend_from_slice(&[2; 64]);
        assert!(matches!(
            File::from_bytes(bytes),
            Err(Hdf5Error::Unsupported(_))
        ));
    }
}
