#ifndef __COGAPS_MUTABLE_MAP_H__
#define __COGAPS_MUTABLE_MAP_H__

#include <stdint.h>

#define ACTIVE_VERSION 1

// VERSION 0: std::map
// VERSION 1: cast away const to modify the key directly in std::map
// VERSION 2: wrap keys in struct with mutable members
// VERSION 3: boost multi-index
// VERSION 4: custom implementation

// This map should provide the same functionality as std::map, with one additional
// feature: it allows the key of an element to be changed, and if that change does
// not re-order the map, then the change should happen in O(1), otherwise the change
// should happen in O(logN) where N is the number of elements in the map

#if ACTIVE_VERSION == 0

#include <map>

template <class K, class V>
class MutableMap
{
public:

    MutableMap() : mMap() {}

    class iterator
    {
    public:
    
        iterator() : mIt() {}
        iterator(const iterator& it) : mIt(it.mIt) {}
        ~iterator() {}
        iterator& operator=(const iterator &it) { mIt = it.mIt; return *this; }
        iterator& operator++() { ++mIt; return *this; }
        iterator& operator--() { --mIt; return *this; }
        std::pair<const K, V>& operator*() { return *mIt; }
        std::pair<const K, V>* operator->() { return &(operator *()); }
        bool operator==(const iterator &it) { return mIt == it.mIt; }
        bool operator!=(const iterator &it) { return !(*this == it); }

    private:

        friend class MutableMap;
        iterator(typename std::map<K, V>::iterator it) : mIt(it) {}
        typename std::map<K, V>::iterator mIt;
    };

    unsigned count(const K &key) const
    {
        return mMap.count(key);
    }
    
    std::pair<iterator, bool> insert(const std::pair<K, V> &val)
    {
        std::pair< typename std::map<K, V>::iterator, bool> result = mMap.insert(val);
        iterator it(result.first);
        return std::pair<iterator, bool>(it, result.second);
    }

    void erase(iterator it)
    {
        mMap.erase(it.mIt);
    }

    iterator find(const K &key)
    {
        return mMap.find(key);
    }

    void updateKey(iterator it, const K &newKey)
    {
        V val = (*it).second;
        mMap.erase(it.mIt);
        insert(std::pair<K, V>(newKey, val));
    }

    iterator begin()
    {
        return iterator(mMap.begin());
    }

    iterator end()
    {
        return iterator(mMap.end());
    }
    
private:

    std::map<K, V> mMap;
};

#elif ACTIVE_VERSION == 1

#include <map>

template <class K, class V>
class MutableMap
{
public:

    MutableMap() : mMap() {}

    class iterator
    {
    public:
    
        iterator() : mIt() {}
        iterator(const iterator& it) : mIt(it.mIt) {}
        ~iterator() {}
        iterator& operator=(const iterator &it) { mIt = it.mIt; return *this; }
        iterator& operator++() { ++mIt; return *this; }
        iterator& operator--() { --mIt; return *this; }
        std::pair<const K, V>& operator*() { return *mIt; }
        std::pair<const K, V>* operator->() { return &(operator *()); }
        bool operator==(const iterator &it) { return mIt == it.mIt; }
        bool operator!=(const iterator &it) { return !(*this == it); }

    private:

        friend class MutableMap;
        iterator(typename std::map<K, V>::iterator it) : mIt(it) {}
        typename std::map<K, V>::iterator mIt;
    };

    class const_iterator
    {
    public:
    
        const_iterator() : mIt() {}
        const_iterator(const const_iterator& it) : mIt(it.mIt) {}
        ~const_iterator() {}
        const_iterator& operator=(const const_iterator &it) { mIt = it.mIt; return *this; }
        const_iterator& operator++() { ++mIt; return *this; }
        const_iterator& operator--() { --mIt; return *this; }
        const std::pair<const K, V>& operator*() { return *mIt; }
        const std::pair<const K, V>* operator->() { return &(operator *()); }
        bool operator==(const const_iterator &it) { return mIt == it.mIt; }
        bool operator!=(const const_iterator &it) { return !(*this == it); }

    private:

        friend class MutableMap;
        const_iterator(typename std::map<K, V>::const_iterator it) : mIt(it) {}
        typename std::map<K, V>::const_iterator mIt;
    };

    unsigned count(const K &key) const
    {
        return mMap.count(key);
    }
    
    std::pair<iterator, bool> insert(const std::pair<K, V> &val)
    {
        std::pair< typename std::map<K, V>::iterator, bool> result = mMap.insert(val);
        iterator it(result.first);
        return std::pair<iterator, bool>(it, result.second);
    }

    void erase(iterator it)
    {
        mMap.erase(it.mIt);
    }

    iterator find(const K &key)
    {
        return mMap.find(key);
    }

    void updateKey(iterator it, const K &newKey)
    {
        const_cast<uint64_t&>((*it).first) = newKey;
    }

    iterator begin()
    {
        return iterator(mMap.begin());
    }

    iterator end()
    {
        return iterator(mMap.end());
    }

    const_iterator begin() const
    {
        return const_iterator(mMap.begin());
    }

    const_iterator end() const
    {
        return const_iterator(mMap.end());
    }
    
private:

    std::map<K, V> mMap;
};

#elif ACTIVE_VERSION == 2

#include <map>

template <class K>
class KeyWrapper
{
public:

    KeyWrapper(K key) : mKey(key) {}

    bool operator<(const KeyWrapper &rhs) const
    {
        return mKey < rhs.mKey;
    }

    void set(K newKey) const
    {
        mKey = newKey;
    }

private:

    mutable K mKey;
};

template <class K, class V>
class MutableMap
{
public:

    MutableMap() : mMap() {}

    class iterator
    {
    public:
    
        iterator() : mIt() {}
        iterator(const iterator& it) : mIt(it.mIt) {}
        ~iterator() {}
        iterator& operator=(const iterator &it) { mIt = it.mIt; return *this; }
        iterator& operator++() { ++mIt; return *this; }
        iterator& operator--() { --mIt; return *this; }
        std::pair<K, V> operator*() { return *mIt; }
        //std::pair<K, V>* operator->() { return &(operator *()); }
        bool operator==(const iterator &it) { return mIt == it.mIt; }
        bool operator!=(const iterator &it) { return !(*this == it); }

    private:

        friend class MutableMap;
        iterator(typename std::map< KeyWrapper<K>, V>::iterator it) : mIt(it) {}
        typename std::map< KeyWrapper<K>, V>::iterator mIt;
    };

    unsigned count(const K &key) const
    {
        return mMap.count(KeyWrapper(key));
    }
    
    std::pair<iterator, bool> insert(const std::pair<K, V> &val)
    {
        std::pair< typename std::map< KeyWrapper<K>, V>::iterator, bool> result = mMap.insert(val);
        iterator it(result.first);
        return std::pair<iterator, bool>(it, result.second);
    }

    void erase(iterator it)
    {
        mMap.erase(it.mIt);
    }

    iterator find(const K &key)
    {
        return mMap.find(KeyWrapper(key));
    }

    void updateKey(iterator it, const K &newKey)
    {
        V val = (*it).second;
        mMap.erase(it.mIt);
        insert(std::pair< KeyWrapper<K>, V>(KeyWrapper(newKey), val));
    }

    iterator begin()
    {
        return iterator(mMap.begin());
    }

    iterator end()
    {
        return iterator(mMap.end());
    }
    
private:

    std::map< KeyWrapper<K>, V> mMap;
};

#endif

#endif // __COGAPS_MUTABLE_MAP_H__

